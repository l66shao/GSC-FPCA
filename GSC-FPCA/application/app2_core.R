# Empirical example, rewritten so that each analysis choice is an option. Same design as the original script unless an
# option changes it. Writes only aggregate summaries and posterior draws of aggregate quantities.
suppressPackageStartupMessages({library(Rcpp); library(RcppArmadillo); library(splines); library(MASS)})
dir.create("cpp_cache", showWarnings = FALSE)
Rcpp::sourceCpp("../method/sampler_tm.cpp", cacheDir = "cpp_cache")
source("../method/loo_cut.R")
source("basis_setup_sparse_Lucy.R")
M <- readRDS("ext_all.rds")
M$sf_bilat <- M$superiorfrontal_grayvol.new      # true left/right average
M$sf_orig <- M$superiorfrontal_grayvol           # column in ncanda_dat.csv; in the version used for the first posting,
                                                 # (left SF volume + left frontal pole vertex count)/2
orig_cols <- c("visit", "sri24_icv_brain_volume", "site", "sex", "visit_age", "race_label", "binge", "id", "base_age", "sexM", "sexF", "superiorfrontal_grayvol", "binge12", "first_binge_visit", "base_age_cat12_14", "base_age_cat15_17", "base_age_cat18_21", "abnormal")
defaults <- list(label = "x", outcome = "sf_bilat", Q = 1, amax = 21, discont = "drop", controls = "all", covs_t = c("icv_baseline", "sexF", "race_labelBlack", "race_labelAsian", "race_labelWhite"), covs_c = NULL,
                 nloop = 2500, burnin = 500, thin = 1, seed = 20261104, psloo = 0, placebo = FALSE, fake = FALSE, keep_draws = TRUE, want = NULL)

prep_data <- function(o) {
  d <- M[complete.cases(M[, orig_cols]), ]
  d <- d[d$abnormal == 0, ]
  d <- d[order(d$id, d$visit), ]
  # later crossing outside the analysis window (for the "no_later" control set)
  crossed_any <- tapply(d$binge12, d$id, max)
  crossed_later <- tapply(d$binge12 * (d$visit_age > o$amax), d$id, max)
  d <- d[d$visit_age >= 12 & d$visit_age <= o$amax, ]
  # outcome must be observed
  d <- d[!is.na(d[[o$outcome]]), ]
  for (v in intersect(c("ses_parent_yoe_b", "fh_alc_density_b", "pds_b"), c(o$covs_t, o$covs_c))) NULL
  ever <- tapply(d$binge12, d$id, max)
  ex_ids <- as.numeric(names(ever)[ever > 0])
  # exposed units whose binge later fell below 12
  mn <- tapply(ifelse(d$binge12 > 0, d$binge, Inf), d$id, min)
  disc_ids <- as.numeric(names(mn)[is.finite(mn) & mn < 12])
  counts <- list(window_ids = length(unique(d$id)), ever = length(ex_ids), disc = length(disc_ids))
  if (o$discont == "drop") {
    d <- d[!d$id %in% disc_ids, ]
  } else if (o$discont == "censor") {
    # keep visits up to the last one before the first report under 12 after onset
    keep <- rep(TRUE, nrow(d))
    for (i in disc_ids) {
      w <- which(d$id == i)
      post <- w[d$binge12[w] > 0]
      lo <- post[d$binge[post] < 12]
      if (length(lo)) keep[w[d$visit[w] >= d$visit[lo[1]]]] <- FALSE
    }
    d <- d[keep, ]
  } # "itt": keep all rows, exposed from first crossing on
  ever <- tapply(d$binge12, d$id, max)
  con_ids <- as.numeric(names(ever)[ever == 0])
  if (o$controls == "no_later") {
    drop <- con_ids[crossed_later[as.character(con_ids)] > 0]
    d <- d[!d$id %in% drop, ]; counts$control_dropped <- length(drop)
  } else if (o$controls == "zero") {
    anyb <- tapply(d$binge > 0, d$id, any)
    drop <- con_ids[anyb[as.character(con_ids)]]
    d <- d[!d$id %in% drop, ]; counts$control_dropped <- length(drop)
  }
  d
}

# Baseline and time-varying covariates built from the first row of each unit in the analysis and from running totals
add_covariates <- function(d, o) {
  d <- d[order(d$id, d$visit), ]
  first <- !duplicated(d$id)
  d$icv_baseline <- ave(d$sri24_icv_brain_volume, d$id, FUN = function(x) x[1])
  d$race_labelWhite <- as.numeric(d$race_label == "Caucasian/White"); d$race_labelAsian <- as.numeric(d$race_label == "Asian")
  d$race_labelBlack <- as.numeric(d$race_label == "African-American/Black"); d$race_labelOther <- as.numeric(d$race_label == "Other")
  for (s in c("B", "C", "D", "E")) d[[paste0("site", s)]] <- as.numeric(d$site == s)
  bl <- function(v) { # value at the visit-1 row of the full data if present, carried to all rows of the unit
    x <- M[[v]][M$visit == 1]; names(x) <- M$id[M$visit == 1]
    as.numeric(x[as.character(d$id)])
  }
  d$ses_yoe_b <- bl("ses_parent_yoe"); d$pds_b <- bl("pds_score"); d$fh_alc_b <- bl("fh_alc_density"); d$fh_drug_b <- bl("fh_drug_density")
  d$mj_life_b <- log1p(ifelse(is.na(bl("youthreport1_cddr55a")), 0, bl("youthreport1_cddr55a")))
  d$tob_life_b <- log1p(ifelse(is.na(bl("youthreport1_cddr4a")), 0, bl("youthreport1_cddr4a")))
  # running lifetime totals as in LifeUseRunningTotal of the extraction script
  rt <- function(a, b) {
    MM <- M[order(M$id, M$visit), ]
    base <- MM[[a]]; inc <- MM[[b]]; inc[is.na(inc)] <- 0; base[is.na(base)] <- 0
    tot <- ave(seq_along(inc), MM$id, FUN = function(k) base[k][1] + cumsum(inc[k]))
    tot[match(paste(d$id, d$visit), paste(MM$id, MM$visit))]
  }
  d$mj_life_t <- log1p(pmax(rt("youthreport1_cddr55a", "youthreport1_cddr55b"), 0))
  d$tob_life_t <- log1p(pmax(rt("youthreport1_cddr4a", "youthreport1_cddr4b"), 0))
  d$zeros <- 0
  d
}

run_app <- function(opt) {
  o <- modifyList(defaults, opt); set.seed(o$seed)
  d <- prep_data(o); d <- add_covariates(d, o)
  if (o$fake) {
    # In-space placebo: drop the exposed, give randomly chosen controls a pseudo-onset drawn from the exposed onset ages.
    # True effect zero, so the estimate should be near zero with an interval that covers it.
    ever <- tapply(d$binge12, d$id, max); ex <- as.numeric(names(ever)[ever > 0])
    ons <- as.numeric(tapply(d$visit_age[d$binge12 > 0], d$id[d$binge12 > 0], min))
    ntr <- length(ex)
    d <- d[!d$id %in% ex, ]
    cand <- sample(as.numeric(names(table(d$id))[table(d$id) >= 2]))
    d$binge12 <- 0; got <- 0
    for (i in cand) {
      if (got >= ntr) break
      w <- which(d$id == i); ag <- d$visit_age[w]
      for (try in 1:20) {
        a0 <- sample(ons, 1)
        k <- which(ag >= a0 & seq_along(ag) >= 2)
        if (length(k)) { d$binge12[w[seq_along(ag) >= k[1]]] <- 1; got <- got + 1; break }
      }
    }
  }
  covs_t <- o$covs_t; covs_c <- if (is.null(o$covs_c)) "zeros" else o$covs_c; est_beta_c <- as.numeric(!is.null(o$covs_c))
  need <- c(covs_t, if (!is.null(o$covs_c)) covs_c)
  miss_units <- unique(d$id[!complete.cases(d[, need, drop = FALSE])])
  d <- d[!d$id %in% miss_units, ]
  # standardize over analysis rows, as in the original script
  for (v in need) if (!v %in% c("sexF", "race_labelBlack", "race_labelAsian", "race_labelWhite", "race_labelOther", "siteB", "siteC", "siteD", "siteE")) d[[v]] <- as.numeric(scale(d[[v]]))
  y0 <- d[[o$outcome]]; center <- mean(y0); sc <- sd(y0); d$y <- (y0 - center) / sc
  d$group <- ifelse(d$binge12 == 1, "exposed", "control")
  if (o$placebo) for (i in unique(d$id[d$binge12 == 1])) {
    w <- which(d$id == i & d$binge12 == 0)
    if (length(w) >= 2) d$group[w[which.max(d$visit_age[w])]] <- "exposed"
  }
  rng <- range(d$visit_age); d$age2 <- (d$visit_age - rng[1]) / (rng[2] - rng[1])
  ctl <- d[d$group == "control", ]; ctl <- ctl[order(ctl$id, ctl$visit), ]
  ids <- sort(unique(ctl$id)); N <- length(ids)
  ex_ids <- sort(unique(d$id[d$group == "exposed"]))
  no_pre <- setdiff(ex_ids, ids)
  bc <- as.numeric(!ids %in% ex_ids)
  sp <- split(ctl, ctl$id)[as.character(ids)]
  Y <- lapply(sp, function(x) matrix(x$y, ncol = 1)); Tm <- lapply(sp, function(x) x$age2)
  Xt <- lapply(sp, function(x) as.matrix(x[, covs_t, drop = FALSE])); Xc <- lapply(sp, function(x) as.matrix(x[, covs_c, drop = FALSE]))
  nknots <- 3; K <- nknots + 4
  bs_ <- basis_setup_sparse(lapply(Tm, list), nknots, delta = 1/1000, plotit = FALSE, orth = FALSE)
  phi_t <- bs_[[3]][[1]]; phi_cont <- bs_[[5]][[1]]
  Q <- o$Q
  if (!is.null(o$holdout)) {
    # Exact leave-one-out for one control: it is marked as treated, so it does not inform the components, and mode 0 returns
    # its log density given the components with prior mean zero for its scores. Nothing else from this fit is used.
    bc2 <- bc; stopifnot(bc[o$holdout] == 1); bc2[o$holdout] <- 0
    f2 <- run_mcmc(o$nloop, o$burnin, o$thin, Y, Tm, Xt, Xc, phi_t, array(0, dim = c(K, 1)), matrix(0, K, length(covs_t)), array(0, length(covs_c)),
                   array(0, dim = c(Q, N)), array(0, dim = c(K, Q)), diag(rep(1, Q), Q), 0.1, diag(rep(1, Q), Q), 1, 1, 1, 1, 1, 1, 1, est_beta_c, 1, bc2, 0)
    return(list(ll = f2[[9]][o$holdout, ]))
  }
  fit <- run_mcmc(o$nloop, o$burnin, o$thin, Y, Tm, Xt, Xc, phi_t, array(0, dim = c(K, 1)), matrix(0, K, length(covs_t)), array(0, length(covs_c)),
                  array(0, dim = c(Q, N)), array(0, dim = c(K, Q)), diag(rep(1, Q), Q), 0.1, diag(rep(1, Q), Q), 1, 1, 1, 1, 1, 1, 1, est_beta_c, o$psloo, bc)
  S <- length(fit[[1]])
  # post-exposure and pre-exposure rows of exposed units
  df <- d[d$id %in% ex_ids & d$id %in% ids, ]; df <- df[order(df$id, df$visit), ]
  lastpre <- tapply(df$visit_age[df$group == "control"], df$id[df$group == "control"], max)
  df$T_scaled <- df$visit_age - lastpre[as.character(df$id)]
  pos <- match(df$id, ids)
  Pm <- phi_cont[, round(df$age2, 3) * 1000 + 1, drop = FALSE]
  Zm <- as.matrix(df[, covs_t, drop = FALSE]); Xm <- as.matrix(df[, covs_c, drop = FALSE])
  Fit0 <- Fit1 <- matrix(NA, nrow(df), S)
  for (s in 1:S) {
    f <- as.vector(crossprod(Pm, fit[[2]][[s]])) + as.vector(Xm %*% fit[[8]][[s]])
    for (j in seq_along(covs_t)) f <- f + Zm[, j] * as.vector(crossprod(Pm, fit[[7]][[s]][, j]))
    TP <- t(Pm) %*% fit[[6]][[s]]
    Fit0[, s] <- f + rowSums(TP * t(fit[[1]][[s]][, pos, drop = FALSE]))
    Fit1[, s] <- f + rowSums(TP * t(fit[[10]][[s]][, pos, drop = FALSE]))
  }
  sg <- unlist(fit[[5]])
  bins <- list()
  if (o$placebo) {
    bins[["placebo"]] <- df$group == "exposed" & df$binge12 == 0
  } else {
    yr <- floor(df$T_scaled + .5); ispost <- df$group == "exposed"
    for (k in 1:5) bins[[sprintf("year %d", k)]] <- ispost & yr == k
    bins[["all"]] <- ispost
    for (k in -3:0) bins[[sprintf("pre %d", k)]] <- !ispost & yr == k
    if (isTRUE(o$age_bins)) {
      bins[["age<=19"]] <- ispost & df$visit_age <= 19; bins[["age 19-20"]] <- ispost & df$visit_age > 19 & df$visit_age <= 20
      bins[["age 20-20.5"]] <- ispost & df$visit_age > 20 & df$visit_age <= 20.5; bins[["age 20.5-20.8"]] <- ispost & df$visit_age > 20.5 & df$visit_age <= 20.8; bins[["age>20.8"]] <- ispost & df$visit_age > 20.8
      bins[["all age<=20.5"]] <- ispost & df$visit_age <= 20.5; bins[["all age<=20.8"]] <- ispost & df$visit_age <= 20.8
    }
  }
  att <- list(); draws <- list()
  for (nm in names(bins)) {
    b <- bins[[nm]]
    if (!any(b)) next
    for (mm in 0:1) {
      Fm <- if (mm == 0) Fit0 else Fit1
      dr <- colMeans(df$y[b] - Fm[b, , drop = FALSE]); nO <- sum(b)
      dp <- dr - rnorm(S, 0, sqrt(sg / nO))
      att[[length(att) + 1]] <- data.frame(label = o$label, method = c("gsc0", "gsc1")[mm + 1], bin = nm, n = nO, units = length(unique(df$id[b])), est = mean(dr), med = median(dr),
        fit_lo = quantile(dr, .025, names = FALSE), fit_hi = quantile(dr, .975, names = FALSE), lo = quantile(dp, .025, names = FALSE), hi = quantile(dp, .975, names = FALSE), post_sd = sd(dr), pred_sd = sd(dp), p_neg = mean(dp < 0), sc = sc)
      draws[[paste(c("gsc0", "gsc1")[mm + 1], nm)]] <- dr
    }
  }
  att <- do.call(rbind, att)
  D_s <- t(sapply(fit[[3]], function(x) diag(as.matrix(x)))); if (Q == 1) D_s <- t(D_s)
  Nu <- t(sapply(fit[[11]], function(x) as.vector(x))); if (Q == 1) Nu <- t(Nu)
  # covariance function draws on an age grid
  ga <- round(seq(0, 1, length.out = 7) * 1000) + 1; Pg <- phi_cont[, ga, drop = FALSE]
  # treated mean expressed on the outcome scale, invariant to rotation and sign of the FPCs: f(t)'nu on the age grid
  NuF <- t(sapply(1:S, function(s) as.vector(t(Pg) %*% fit[[6]][[s]] %*% fit[[11]][[s]])))
  Cv <- t(sapply(1:S, function(s) { Ff <- t(Pg) %*% fit[[6]][[s]]; Cm <- Ff %*% fit[[3]][[s]] %*% t(Ff); c(diag(Cm), Cm[1, 7], Cm[1, 4], Cm[4, 7]) }))
  gi <- seq(1, 1001, by = 5); Pc <- phi_cont[, gi, drop = FALSE]; Cbar <- matrix(0, length(gi), length(gi))
  for (s_ in 1:S) { Ff <- t(Pc) %*% fit[[6]][[s_]]; Cbar <- Cbar + Ff %*% fit[[3]][[s_]] %*% t(Ff) / S }
  ec <- eigen(Cbar, symmetric = TRUE)
  meta <- data.frame(label = o$label, outcome = o$outcome, Q = Q, amax = o$amax, discont = o$discont, controls = o$controls, covs_t = paste(covs_t, collapse = "+"), covs_c = paste(if (is.null(o$covs_c)) "" else covs_c, collapse = "+"),
    N = N, Ntr = sum(bc == 0), Nco = sum(bc == 1), no_pre = length(no_pre), miss_units = length(miss_units), rows = nrow(d), rows_ctl = nrow(ctl), n_post = sum(df$group == "exposed"), center = center, sc = sc, rng_lo = rng[1], rng_hi = rng[2],
    sigma2 = mean(sg), sigma2_lo = quantile(sg, .025, names = FALSE), sigma2_hi = quantile(sg, .975, names = FALSE), d1 = mean(D_s[, 1]), nu1 = mean(Nu[, 1]), nu1_lo = quantile(Nu[, 1], .025, names = FALSE), nu1_hi = quantile(Nu[, 1], .975, names = FALSE),
    nu1_sd = mean(Nu[, 1] / sqrt(D_s[, 1])), nu1_sd_lo = quantile(Nu[, 1] / sqrt(D_s[, 1]), .025, names = FALSE), nu1_sd_hi = quantile(Nu[, 1] / sqrt(D_s[, 1]), .975, names = FALSE), S = S, seed = o$seed, placebo = o$placebo)
  out <- list(att = att, meta = meta, draws = draws, sigma2 = sg, D = D_s, Nu = Nu, Cov = Cv, NuF = NuF, age_grid = rng[1] + (rng[2] - rng[1]) * (ga - 1) / 1000, eig_share = ec$values[1:Q] / sum(ec$values), eig_fun = ec$vectors[, 1:Q, drop = FALSE] / sqrt(5 / 1000), eig_age = rng[1] + (rng[2] - rng[1]) * (gi - 1) / 1000)
  if (o$psloo == 1) { out$loo <- elpd_cut(fit[[9]], bc == 1) }
  saveRDS(out, file.path("res2", paste0(o$label, ".rds")))
  write.csv(att, file.path("res2", paste0(o$label, "_att.csv")), row.names = FALSE)
  out
}
