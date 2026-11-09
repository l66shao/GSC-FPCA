# Known components: coverage of predictive intervals when treated scores and errors are not normal.
# Checks the Lyapunov results (asymptotic normality of the sampling law) for prior mean zero and a separate prior mean.
set.seed(20261106)
f <- function(t) cbind(sqrt(2) * sin(2 * pi * t), sqrt(2) * cos(2 * pi * t))
D0 <- diag(c(2, 1)); s2 <- 1; Dinv <- solve(D0); k <- 2; z <- qnorm(.975)
make_design <- function(Ntr, m) lapply(1:Ntr, function(i) {
  tp <- sort(runif(m, 0, .6)); npost <- sample(1:2, 1); tq <- sort(runif(npost, .6, 1))
  Fi <- f(tp); ci <- colSums(f(tq))
  Prec <- crossprod(Fi) / s2 + Dinv; Om <- solve(Prec); R <- Om %*% Dinv; W <- Om %*% t(Fi) / s2
  list(F = Fi, c = ci, n = npost, Om = Om, R = R, W = W, m = m)
})
# Selected scores: population N(0, D0), treated if alpha1 + alpha2 + logistic noise > 1.5.
thr <- 1.5; LD <- chol(D0)
sel_draw <- function(n) {
  out <- NULL; got <- 0
  while (got < n) { a <- matrix(rnorm(2 * 5 * n), ncol = 2) %*% LD
    a <- a[a[, 1] + a[, 2] + rlogis(nrow(a)) > thr, , drop = FALSE]; out <- rbind(out, a); got <- nrow(out) }
  out[1:n, , drop = FALSE]
}
big <- sel_draw(4e6); nu_sel <- colMeans(big); Sig_sel <- cov(big)
psel <- mean(local({a <- matrix(rnorm(2e6), ncol = 2) %*% LD; a[, 1] + a[, 2] + rlogis(nrow(a)) > thr}))
rm(big)
std <- function(n, type) switch(type, normal = rnorm(n), skewed = rexp(n) - 1, t3 = rt(n, 3) / sqrt(3))
scores <- function(nrep, type, nu, Sig) {
  if (type == "selected") return(sel_draw(nrep))
  matrix(nu, nrep, k, byrow = TRUE) + matrix(std(nrep * k, type), nrep, k) %*% chol(Sig)
}
errs_fun <- function(n, type) sqrt(s2) * std(n, if (type == "selected") "normal" else type)
run <- function(Ntr, m, type, nrep = 20000) {
  des <- make_design(Ntr, m); nO <- sum(sapply(des, `[[`, "n"))
  nu <- if (type == "selected") nu_sel else c(1, 1); Sig <- if (type == "selected") Sig_sel else D0
  G <- Reduce(`+`, lapply(des, function(x) Dinv %*% (diag(k) - x$R))); Ginv <- solve(G)
  ell <- Reduce(`+`, lapply(des, function(x) t(x$R) %*% x$c)) / nO
  vc0 <- sum(sapply(des, function(x) t(x$c) %*% x$Om %*% x$c)) / nO^2
  v0 <- vc0 + s2 / nO; vnu <- v0 + as.numeric(t(ell) %*% Ginv %*% ell); vcnu <- vnu - s2 / nO
  B0 <- sum(sapply(des, function(x) t(x$c) %*% x$R %*% nu)) / nO; Bnu <- 0
  y <- Dinv %*% Ginv %*% ell
  vcS0 <- sum(sapply(des, function(x) t(x$c) %*% (s2 * x$W %*% t(x$W) + x$R %*% Sig %*% t(x$R)) %*% x$c)) / nO^2
  vS0 <- vcS0 + s2 / nO
  vcSnu <- sum(sapply(des, function(x) { a <- t(x$R) %*% x$c / nO - t(diag(k) - x$R) %*% y; b <- -t(x$W) %*% (x$c / nO + y)
    as.numeric(t(a) %*% Sig %*% a + s2 * t(b) %*% b) }))
  vSnu <- vcSnu + s2 / nO
  bsum <- matrix(0, nrep, k); e0 <- numeric(nrep); epsbar <- numeric(nrep)
  for (x in des) {
    al <- scores(nrep, type, nu, Sig)
    u <- al %*% t(x$F) + matrix(errs_fun(nrep * x$m, type), nrep, x$m)
    a <- u %*% t(x$W)
    bsum <- bsum + a %*% Dinv
    ep <- rowSums(matrix(errs_fun(nrep * x$n, type), nrep, x$n))
    e0 <- e0 + as.vector((al - a) %*% x$c) + ep; epsbar <- epsbar + ep
  }
  nuh <- bsum %*% t(Ginv)
  e1 <- e0 - as.vector(nuh %*% (ell * nO))
  e0 <- e0 / nO; e1 <- e1 / nO; epsbar <- epsbar / nO
  pif <- function(B, sv, sS) pnorm((z * sv - B) / sS) - pnorm((-z * sv - B) / sS)
  data.frame(type = type, Ntr = Ntr, m = m, nO = nO, B0 = B0,
    cov0_formula = pif(B0, sqrt(v0), sqrt(vS0)), cov0_mc = mean(abs(e0) <= z * sqrt(v0)),
    covnu_formula = pif(Bnu, sqrt(vnu), sqrt(vSnu)), covnu_mc = mean(abs(e1) <= z * sqrt(vnu)),
    covnu_fit_formula = pif(Bnu, sqrt(vcnu), sqrt(vcSnu)), covnu_fit_mc = mean(abs(e1 - epsbar) <= z * sqrt(vcnu)),
    sd_formula = sqrt(vSnu), sd_mc = sd(e1), skew_mc = mean((e1 - mean(e1))^3) / sd(e1)^3,
    exkurt_mc = mean((e1 - mean(e1))^4) / sd(e1)^4 - 3)
}
res <- list()
for (type in c("normal", "skewed", "t3", "selected")) for (m in c(1, 4)) for (Ntr in c(20, 80, 320))
  res[[length(res) + 1]] <- run(Ntr, m, type)
out <- do.call(rbind, res)
write.csv(out, "nongauss_check.csv", row.names = FALSE)
cat(sprintf("selection prob %.3f; nu_sel = (%.3f, %.3f); Sig_sel = [%.3f %.3f; %.3f %.3f]\n", psel, nu_sel[1], nu_sel[2], Sig_sel[1,1], Sig_sel[1,2], Sig_sel[2,1], Sig_sel[2,2]))
print(cbind(out[, 1:4], round(out[, -(1:4)], 3)), row.names = FALSE)
