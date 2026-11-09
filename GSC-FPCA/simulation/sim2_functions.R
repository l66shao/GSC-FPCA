# GSC-FPCA expanded simulation: data-generating scenarios and competing estimators.
# Sourced by run_sim2.R. See the header of run_sim2.R for the design.
library(splines); library(MASS)

grid20 <- seq(0, 1, 1/19)

# spline basis on the 20-point grid (cubic B-splines, 3 interior knots, each column scaled as in the fit)
truth_basis <- function(q_inner = 3) {
  b <- basis_setup_sparse(list(list(grid20)), q_inner, delta = 1/600)
  t(b[[5]][[1]][, round(grid20 * 600) + 1, drop = FALSE])  # 20 x q
}

# True curves on the 20-point grid for each scenario. Returns mu (20), Fm (20 x k), Gam (20) and score sds.
truth_curves <- function(dgp, params) {
  B <- truth_basis()
  mu <- as.vector(B %*% params[[7]])
  Gam <- as.vector(B %*% rep(1, 7))
  Fm <- B %*% params[[6]]
  if (dgp == "outside") {            # curves outside the span of the fitted basis
    mu <- 2 * sin(3 * pi * grid20) + 2 * cos(5 * pi * grid20)
    Fm <- cbind(sqrt(2) * sin(4 * pi * grid20), sqrt(2) * cos(4 * pi * grid20))
  } else if (dgp == "wrongk") {      # three components, fitted with two
    x <- seq(0, 2 * pi, length.out = 7)
    th3 <- sin(2 * x); th3 <- th3 / sqrt(sum(th3^2))
    Fm <- cbind(Fm, B %*% th3)
  } else if (dgp == "growth") {      # intercept and slope, not generated from FPCs
    Fm <- cbind(rep(1, 20), 2 * (grid20 - .5))
  }
  list(mu = mu, Fm = Fm, Gam = Gam)
}

# Simulate one data set. Returns a long data frame with one row per visit plus unit-level fields.
# Selection on scores: U = b'alpha + logistic noise, and the N_tr largest U are treated, so N_tr is fixed.
simulate2 <- function(Ntr, Nco, lambda, n_min, dgp, tc, select_b = 1, n_cap = 10) {
  N <- Ntr + Nco
  k <- ncol(tc$Fm)
  # scores
  if (dgp == "outside") {            # skewed scores: centered and scaled Gamma(2)
    A <- matrix((rgamma(N * k, 2) - 2) / sqrt(2), N, k)
  } else if (dgp == "growth") {      # intercept and slope, correlation 0.3
    A <- mvrnorm(N, c(0, 0), matrix(c(1, .3, .3, 1), 2))
  } else {
    A <- matrix(rnorm(N * k), N, k)
  }
  # who is treated
  if (dgp == "none") {
    trt <- c(rep(FALSE, Nco), rep(TRUE, Ntr))
  } else {
    U <- as.vector(A %*% rep(select_b, k)) + rlogis(N)
    trt <- rank(-U, ties.method = "first") <= Ntr
  }
  rows <- vector("list", N)
  for (i in 1:N) {
    lam_i <- lambda
    if (dgp == "visits") lam_i <- lambda * exp(.4 * A[i, 1])  # visit counts depend on the first score
    n <- min(max(rpois(1, lam_i), n_min), n_cap)
    idx <- sort(sample(20)[1:n])
    t_i <- grid20[idx]
    Z1 <- rnorm(1); Z2 <- rbinom(1, 1, .5)
    X1 <- rnorm(n); X2 <- rbinom(n, 1, .5)
    # errors
    if (dgp == "outside") {
      e <- rt(n, 5) * sqrt(3 / 5)     # t_5 with variance 1
    } else if (dgp == "growth") {     # Ornstein-Uhlenbeck plus white noise, total variance 1
      R <- exp(-abs(outer(t_i, t_i, "-")) / .2)
      e <- as.vector(t(chol(R + 1e-10 * diag(n))) %*% rnorm(n)) * sqrt(.5) + rnorm(n, 0, sqrt(.5))
    } else {
      e <- rnorm(n)
    }
    y0 <- tc$mu[idx] + (Z1 + Z2) * tc$Gam[idx] + X1 + X2 + as.vector(tc$Fm[idx, , drop = FALSE] %*% A[i, ]) + e
    m <- n
    if (trt[i]) {
      m <- max(floor(2 * n / 3), 1)
      if (dgp == "dip" && n >= 2) {   # exposure follows a transient low: T0 is the visit with the lowest error
        m <- which.min(e[1:(n - 1)])
      }
    }
    delta <- ifelse(seq_len(n) > m, t_i - t_i[m], NA_real_)
    rows[[i]] <- data.frame(id = i, time_idx = idx, t = t_i, y0 = y0, e = e, X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2,
                            trt = trt[i], post = trt[i] & seq_len(n) > m, delta = delta, visit = seq_len(n))
  }
  d <- do.call(rbind, rows)
  attr(d, "nu") <- colMeans(A[trt, , drop = FALSE]) - colMeans(A[!trt, , drop = FALSE])
  d
}

effect_fun <- function(delta, shape) {
  switch(shape,
    linear = 5 * delta,
    constant = rep(1, length(delta)),
    saturating = 2 * (1 - exp(-5 * delta)),
    transient = 10 * delta * exp(-5 * delta))
}

# Fit GSC-FPCA once. The chain gives both prior-mean versions. fit[[12]] holds the log likelihood when psloo = 1.
fit_gsc <- function(d, Q, nloop = 1000, burnin = 500, psloo = 0) {
  ids <- sort(unique(d$id)); N <- length(ids)
  Y <- T_ <- X <- Xc <- vector("list", N); bc <- numeric(N)
  for (i in ids) {
    di <- d[d$id == i & !d$post, ]
    Y[[i]] <- list(matrix(di$y0)); T_[[i]] <- list(matrix(di$t))
    x <- cbind(X1 = di$Z1, X2 = di$Z2); xc <- cbind(X1 = di$X1, X2 = di$X2)
    X[[i]] <- list(x); Xc[[i]] <- list(xc); bc[i] <- as.numeric(!di$trt[1])
  }
  mfpca_1d(Y_sparse = Y, time_sparse = T_, X_sparse = X, X_c_sparse = Xc, group = rep(1, N), Q = Q, nknots = 3,
           nloop = nloop, burnin = burnin, thin = 1, orth = FALSE, sim = FALSE, PSLOO = psloo, b_control = bc)
}

# Membership of post-exposure observations in a reporting bin. Bin labels are event time rounded to 0.1, or "all".
in_bin <- function(b, binvec) if (b == "all") rep(TRUE, length(binvec)) else binvec == b

# ATT errors for GSC-FPCA under both prior means. An error is observed no-treatment outcome minus counterfactual,
# averaged over the bin; adding the true effect gives the estimate, so the error of every counterfactual-based
# estimator is the same for any effect shape. Intervals are for the error and cover the truth when they contain 0.
gsc_att <- function(fit, d, bins) {
  po <- d[d$post, ]
  pc <- fit[[9]][[5]][[1]]
  A0 <- fit[[1]]; A1 <- fit[[10]]; Mu <- fit[[2]]; Sg <- unlist(fit[[5]]); Th <- fit[[6]]; Be <- fit[[7]]; Bc <- fit[[8]]; Nu <- fit[[11]]
  S <- length(A0)
  Pm <- pc[, round(po$t, 3) * 1000 + 1, drop = FALSE]
  Xm <- cbind(po$Z1, po$Z2); Xcm <- cbind(po$X1, po$X2)
  Y0 <- matrix(NA, nrow(po), S); Y1 <- Y0
  for (s in 1:S) {
    f <- crossprod(Pm, Mu[[s]]) + Xcm %*% Bc[[s]]
    for (j in 1:2) f <- f + Xm[, j] * crossprod(Pm, Be[[s]][, j])
    TP <- t(Pm) %*% Th[[s]]
    Y0[, s] <- f + rowSums(TP * t(A0[[s]][, po$id, drop = FALSE]))
    Y1[, s] <- f + rowSums(TP * t(A1[[s]][, po$id, drop = FALSE]))
  }
  Nz <- matrix(rnorm(length(Y0)), nrow(Y0), S) * rep(sqrt(Sg), each = nrow(Y0))
  out <- list()
  for (b in bins) {
    w <- in_bin(b, po$bin)
    if (!any(w)) next
    for (mm in 0:1) {
      Ym <- if (mm == 0) Y0 else Y1
      draws <- colMeans(po$y0[w] - Ym[w, , drop = FALSE])
      a <- draws - colMeans(Nz[w, , drop = FALSE])
      out[[length(out) + 1]] <- data.frame(method = if (mm == 0) "gsc0" else "gsc1", bin = b, n = sum(w),
        truth_linear = mean(5 * po$delta[w]), err = mean(draws), err_med = median(draws),
        lo = quantile(a, .025, names = FALSE), hi = quantile(a, .975, names = FALSE),
        indiv_mse = mean((po$y0[w] - po$e[w] - rowMeans(Ym[w, , drop = FALSE]))^2))
    }
  }
  r <- do.call(rbind, out)
  r$sigma2 <- mean(Sg)
  r
}

# Difference-in-differences by imputation: unit effects plus a common time curve and covariates, fit on untreated
# observations (Borusyak, Jaravel and Spiess), with a cluster bootstrap over units. Same spline for the time curve.
did_att <- function(d, B_mat, bins, nboot = 200) {
  Wb <- B_mat[d$time_idx, , drop = FALSE]
  W <- cbind(Wb, d$Z1 * Wb, d$Z2 * Wb, d$X1, d$X2)
  ids <- sort(unique(d$id)); p <- ncol(W); n_u <- length(ids)
  AAv <- matrix(0, p * p, n_u); cc <- matrix(0, p, n_u); res <- vector("list", n_u)
  groups <- split(seq_len(nrow(d)), d$id)
  for (i in seq_along(ids)) {
    r <- groups[[as.character(ids[i])]]
    ru <- r[!d$post[r]]
    wbar <- colMeans(W[ru, , drop = FALSE]); ybar <- mean(d$y0[ru])
    Wt <- sweep(W[ru, , drop = FALSE], 2, wbar); yt <- d$y0[ru] - ybar
    AAv[, i] <- as.vector(crossprod(Wt)); cc[, i] <- crossprod(Wt, yt)
    rp <- r[d$post[r]]
    if (length(rp)) res[[i]] <- list(rp = rp, y = d$y0[rp] - ybar, q = sweep(W[rp, , drop = FALSE], 2, wbar), e = d$e[rp])
  }
  ridge <- 1e-8 * diag(p)
  th_fun <- function(om) solve(matrix(AAv %*% om, p, p) + ridge, as.vector(cc %*% om))
  tr <- which(!vapply(res, is.null, TRUE))
  nb <- length(bins)
  P <- matrix(0, length(tr), nb); Nn <- P; Q <- lapply(1:nb, function(b) matrix(0, length(tr), p))
  for (u in seq_along(tr)) {
    x <- res[[tr[u]]]
    for (b in 1:nb) {
      w <- in_bin(bins[b], d$bin[x$rp])
      if (any(w)) { P[u, b] <- sum(x$y[w]); Nn[u, b] <- sum(w); Q[[b]][u, ] <- colSums(x$q[w, , drop = FALSE]) }
    }
  }
  est_fun <- function(th_, om_tr) vapply(1:nb, function(b) sum(om_tr * (P[, b] - Q[[b]] %*% th_)) / sum(om_tr * Nn[, b]), 0)
  th_hat <- th_fun(rep(1, n_u))
  est <- est_fun(th_hat, rep(1, length(tr)))
  # squared error of each unit's counterfactual mean: y - e - (ybar + q th)
  Sq <- numeric(nb)
  for (u in seq_along(tr)) {
    x <- res[[tr[u]]]
    sq <- (x$y - x$e - as.vector(x$q %*% th_hat))^2
    for (b in 1:nb) { w <- in_bin(bins[b], d$bin[x$rp]); if (any(w)) Sq[b] <- Sq[b] + sum(sq[w]) }
  }
  bs <- matrix(NA, nboot, nb)
  for (r in 1:nboot) {
    om <- tabulate(sample.int(n_u, replace = TRUE), n_u)
    bs[r, ] <- est_fun(th_fun(om), om[tr])
  }
  ok <- colSums(Nn) > 0
  r <- data.frame(method = "did", bin = bins, n = colSums(Nn), truth_linear = NA_real_, err = est, err_med = NA_real_,
             lo = apply(bs, 2, quantile, .025, na.rm = TRUE), hi = apply(bs, 2, quantile, .975, na.rm = TRUE),
             indiv_mse = Sq / pmax(colSums(Nn), 1), sigma2 = NA_real_)
  r[ok, ]
}

# Linear mixed models, the common practice: spline time curve, covariates, random intercept and slope.
# lme1 has one exposure coefficient, as in the empirical example. lme2 has one coefficient per 0.1 event-time bin,
# and its "all" estimate is their observation-weighted mean. Returns errors relative to the truth for each effect shape.
lme_att <- function(d, B_mat, bins, shapes = "linear") {
  W <- B_mat[d$time_idx, , drop = FALSE]; colnames(W) <- paste0("b", 1:ncol(W))
  ZW1 <- d$Z1 * W; ZW2 <- d$Z2 * W; colnames(ZW1) <- paste0("z1", 1:ncol(W)); colnames(ZW2) <- paste0("z2", 1:ncol(W))
  X <- cbind(W[, -1], ZW1, ZW2, X1 = d$X1, X2 = d$X2)
  ctrl <- lme4::lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore", check.conv.grad = "ignore", check.conv.hess = "ignore",
                      check.rankX = "silent.drop.cols", check.scaleX = "ignore", optimizer = "bobyqa")
  out <- list()
  lv <- sort(unique(d$bin[d$post]))
  Pd <- sapply(lv, function(b) as.numeric(d$post & !is.na(d$bin) & d$bin == b)); colnames(Pd) <- paste0("p", gsub("[.]", "", lv))
  for (shape in shapes) {
    y <- d$y0
    y[d$post] <- y[d$post] + effect_fun(d$delta[d$post], shape)
    tr_all <- mean(effect_fun(d$delta[d$post], shape))
    dd <- data.frame(y = y, id = d$id, t = d$t, post = as.numeric(d$post), X)
    m1 <- tryCatch(suppressWarnings(suppressMessages(lme4::lmer(reformulate(c(colnames(X), "post", "(1 + t | id)"), response = "y"), data = dd, REML = TRUE, control = ctrl))), error = function(e) NULL)
    if (!is.null(m1)) {
      cs <- coef(summary(m1))["post", ]
      out[[length(out) + 1]] <- data.frame(method = "lme1", shape = shape, bin = "all", n = sum(d$post), err = cs[1] - tr_all, lo = cs[1] - tr_all - 1.96 * cs[2], hi = cs[1] - tr_all + 1.96 * cs[2])
    }
    dd2 <- data.frame(dd[, setdiff(names(dd), "post")], Pd)
    m2 <- tryCatch(suppressWarnings(suppressMessages(lme4::lmer(reformulate(c(colnames(X), colnames(Pd), "(1 + t | id)"), response = "y"), data = dd2, REML = TRUE, control = ctrl))), error = function(e) NULL)
    if (!is.null(m2)) {
      fe <- lme4::fixef(m2); nmz <- intersect(colnames(Pd), names(fe))
      be <- fe[nmz]; V <- as.matrix(vcov(m2))[nmz, nmz, drop = FALSE]
      nb <- colSums(Pd)[nmz]; wt <- nb / sum(nb)
      allest <- sum(wt * be); allse <- sqrt(as.numeric(t(wt) %*% V %*% wt))
      out[[length(out) + 1]] <- data.frame(method = "lme2", shape = shape, bin = "all", n = sum(d$post), err = allest - tr_all, lo = allest - tr_all - 1.96 * allse, hi = allest - tr_all + 1.96 * allse)
      for (b in setdiff(bins, "all")) {
        nm <- paste0("p", gsub("[.]", "", b))
        if (nm %in% nmz) {
          w <- d$post & !is.na(d$bin) & d$bin == b
          trb <- mean(effect_fun(d$delta[w], shape))
          e <- be[[nm]] - trb; se <- sqrt(V[nm, nm])
          out[[length(out) + 1]] <- data.frame(method = "lme2", shape = shape, bin = b, n = sum(w), err = e, lo = e - 1.96 * se, hi = e + 1.96 * se)
        }
      }
    }
  }
  r <- do.call(rbind, out); rownames(r) <- NULL
  r
}

# Interactive fixed effects (fect, method "ife") and matrix completion (method "mc") on the 20-point grid, with the
# same covariates as the other estimators. The fit imputes from untreated cells; the simulated Y is the no-treatment
# outcome, so eff on a post-exposure cell is the error of the estimate there. Units with fewer than min_T0 untreated
# cells are dropped by the package; errors are averaged over the cells it keeps, and "kept" records the share of
# treated units retained. r is fixed unless cv = TRUE. With se = TRUE the bootstrap interval is for "all" only.
fect_att <- function(d, B_mat, bins, method = "ife", r = 2, min_T0 = NULL, cv = NULL, rmax = 3, se = FALSE, nboots = 100, label = method) {
  W <- B_mat[d$time_idx, ]
  Zc <- cbind(d$Z1 * W[, -1], d$Z2 * W[, -1]); colnames(Zc) <- paste0("zb", 1:ncol(Zc))
  dd <- data.frame(Y = d$y0, D = as.numeric(d$post), id = d$id, time = d$time_idx, X1 = d$X1, X2 = d$X2, Zc)
  covs <- c("X1", "X2", colnames(Zc))
  fm <- as.formula(paste("Y ~ D +", paste(covs, collapse = "+")))
  if (is.null(cv)) cv <- (method == "mc")
  args <- list(formula = fm, data = dd, index = c("id", "time"), force = "two-way", method = method, CV = cv, se = se, parallel = FALSE, na.rm = FALSE)
  if (method == "ife") args$r <- if (cv) c(0, rmax) else r
  if (cv) { args$k <- 3 }
  if (method == "mc") args$nlambda <- 5
  if (!is.null(min_T0)) args$min.T0 <- min_T0
  if (se) { args$nboots <- nboots; args$vartype <- "bootstrap" }
  invisible(capture.output(fit <- tryCatch(suppressWarnings(suppressMessages(do.call(fect::fect, args))), error = function(e) NULL)))
  if (is.null(fit)) return(NULL)
  po <- d[d$post, ]
  i1 <- match(po$time_idx, fit$rawtime); i2 <- match(po$id, fit$id); ok <- !is.na(i2)
  if (!any(ok)) return(NULL)
  e <- fit$eff[cbind(i1[ok], i2[ok])]
  yct <- fit$Y.ct[cbind(i1[ok], i2[ok])]
  po <- po[ok, ]
  kept <- length(unique(po$id)) / length(unique(d$id[d$trt]))
  out <- list()
  for (b in bins) {
    w <- in_bin(b, po$bin)
    if (!any(w)) next
    out[[length(out) + 1]] <- data.frame(method = label, bin = b, n = sum(w), truth_linear = mean(5 * po$delta[w]), err = mean(e[w]), err_med = NA_real_,
      lo = NA_real_, hi = NA_real_, indiv_mse = mean((po$y0[w] - po$e[w] - yct[w])^2), sigma2 = NA_real_)
  }
  r_out <- do.call(rbind, out)
  r_out$kept <- kept
  r_out$r_used <- if (method == "ife") fit$r.cv else NA_real_
  if (se) {
    a <- fit$est.avg
    j <- which(r_out$bin == "all")
    r_out$lo[j] <- a[1, "CI.lower"]; r_out$hi[j] <- a[1, "CI.upper"]
  }
  r_out
}
