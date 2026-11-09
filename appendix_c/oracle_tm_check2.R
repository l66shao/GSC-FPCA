# Known components: Monte Carlo check of Theorem (exact coverage, separate prior mean) including the general variance v_Sigma.
set.seed(20261104)
f <- function(t) cbind(sqrt(2) * sin(2 * pi * t), sqrt(2) * cos(2 * pi * t))
D0 <- diag(c(2, 1)); s2 <- 1; Dinv <- solve(D0); k <- 2
make_design <- function(Ntr, m) lapply(1:Ntr, function(i) {
  tp <- sort(runif(m, 0, .6)); npost <- sample(1:2, 1); tq <- sort(runif(npost, .6, 1))
  Fi <- f(tp); ci <- colSums(f(tq))
  Prec <- crossprod(Fi) / s2 + Dinv; Om <- solve(Prec); R <- Om %*% Dinv; W <- Om %*% t(Fi) / s2
  list(F = Fi, c = ci, n = npost, Om = Om, R = R, W = W, m = m)
})
run <- function(Ntr, m, nu_fun, Sig = D0, nrep = 4000) {
  des <- make_design(Ntr, m); nO <- sum(sapply(des, `[[`, "n"))
  G <- Reduce(`+`, lapply(des, function(x) Dinv %*% (diag(k) - x$R))); Ginv <- solve(G)
  ell <- Reduce(`+`, lapply(des, function(x) t(x$R) %*% x$c)) / nO
  v0 <- sum(sapply(des, function(x) t(x$c) %*% x$Om %*% x$c)) / nO^2 + s2 / nO
  vnu <- v0 + as.numeric(t(ell) %*% Ginv %*% ell)
  nus <- lapply(seq_along(des), nu_fun)
  nustar <- Ginv %*% Reduce(`+`, Map(function(x, nu) Dinv %*% (diag(k) - x$R) %*% nu, des, nus))
  Bnu <- sum(unlist(Map(function(x, nu) t(x$c) %*% x$R %*% (nu - nustar), des, nus))) / nO
  B0 <- sum(unlist(Map(function(x, nu) t(x$c) %*% x$R %*% nu, des, nus))) / nO
  y <- Dinv %*% Ginv %*% ell
  vS <- s2 / nO + sum(sapply(des, function(x) { a <- t(x$R) %*% x$c / nO - t(diag(k) - x$R) %*% y; b <- -t(x$W) %*% (x$c / nO + y); as.numeric(t(a) %*% Sig %*% a + s2 * t(b) %*% b) }))
  LS <- chol(Sig); errs <- matrix(NA, nrep, 2)
  for (r in 1:nrep) {
    b <- 0; a <- vector("list", Ntr); eps_post <- 0; alpha <- vector("list", Ntr)
    for (i in 1:Ntr) { x <- des[[i]]
      alpha[[i]] <- nus[[i]] + as.vector(t(LS) %*% rnorm(k))
      u <- x$F %*% alpha[[i]] + rnorm(x$m, 0, sqrt(s2))
      a[[i]] <- as.vector(x$W %*% u); b <- b + Dinv %*% a[[i]]
      eps_post <- eps_post + sum(rnorm(x$n, 0, sqrt(s2))) }
    nuh <- as.vector(Ginv %*% b); e0 <- 0; e1 <- 0
    for (i in 1:Ntr) { x <- des[[i]]
      e0 <- e0 + sum(x$c * (alpha[[i]] - a[[i]]))
      e1 <- e1 + sum(x$c * (alpha[[i]] - a[[i]] - as.vector(x$R %*% nuh))) }
    errs[r, ] <- (c(e0, e1) + eps_post) / nO
  }
  z <- qnorm(.975); sv <- sqrt(vnu); sS <- sqrt(vS)
  data.frame(Ntr = Ntr, m = m, nO = nO, B0_formula = B0, B0_mc = mean(errs[, 1]), Bnu_formula = Bnu, Bnu_mc = mean(errs[, 2]), se_mc = sd(errs[, 2]) / sqrt(nrep),
             sd_post = sv, sd_formula = sS, sd_mc = sd(errs[, 2]),
             cov_formula = pnorm((z * sv - Bnu) / sS) - pnorm((-z * sv - Bnu) / sS), cov_mc = mean(abs(errs[, 2]) <= z * sv))
}
res <- list()
for (m in c(1, 2, 4, 8)) for (Ntr in c(40, 160)) res[[length(res) + 1]] <- cbind(case = "common nu=(1,1)", run(Ntr, m, function(i) c(1, 1)))
for (m in c(1, 4)) res[[length(res) + 1]] <- cbind(case = "nu=(0,0)", run(40, m, function(i) c(0, 0)))
for (m in c(1, 4)) res[[length(res) + 1]] <- cbind(case = "mixed nu", run(80, m, function(i) if (i %% 2) c(2, 0) else c(0, 2)))
for (m in c(1, 4)) res[[length(res) + 1]] <- cbind(case = "nu=(1,1), Sigma=2*D0", run(80, m, function(i) c(1, 1), Sig = 2 * D0))
for (m in c(1, 4)) res[[length(res) + 1]] <- cbind(case = "nu=(1,1), Sigma=D0/2", run(80, m, function(i) c(1, 1), Sig = D0 / 2))
out <- do.call(rbind, res)
write.csv(out, "oracle_tm_check2.csv", row.names = FALSE)
print(cbind(out[, 1:4], round(out[, -(1:4)], 3)), row.names = FALSE)
