# Compare the sampler's leave-one-out densities with an independent calculation from stored draws, and check that the chain is unchanged.
library(Rcpp); library(MASS); library(splines); pdf(NULL); Sys.setenv(COMEGA = "1")
Rcpp::sourceCpp("../method/sampler_tm.cpp", cacheDir = "cpp_cache")
source("../method/mfpca.R"); source("base_functions.R"); source("sim2_functions.R"); source("../method/loo_cut.R")
params <- make_params(); set.seed(5)
tc <- truth_curves("select", params); d <- simulate2(40, 100, 4, 3, "select", tc); d$bin <- ifelse(d$post, sprintf("%.1f", round(d$delta * 10) / 10), NA)
set.seed(11); f <- fit_gsc(d, 2, nloop = 300, burnin = 200, psloo = 1)
ll <- f[[12]]; dim(ll); ids <- sort(unique(d$id)); trt <- sapply(ids, function(i) d$trt[d$id == i][1])
pc <- f[[9]][[5]][[1]]; s <- 77
Mu <- f[[2]][[s]]; Th <- f[[6]][[s]]; Be <- f[[7]][[s]]; Bc <- f[[8]][[s]]; Dm <- f[[3]][[s]]; sg <- as.numeric(f[[5]][[s]])
pieces <- lapply(which(trt), function(i) { di <- d[d$id == i & !d$post, ]; P <- pc[, round(di$t, 3) * 1000 + 1, drop = FALSE]
  g <- crossprod(P, Mu) + cbind(di$X1, di$X2) %*% Bc; g <- g + di$Z1 * crossprod(P, Be[, 1]) + di$Z2 * crossprod(P, Be[, 2])
  u <- di$y0 - as.vector(g); Fi <- t(P) %*% Th; V <- Fi %*% Dm %*% t(Fi) + sg * diag(length(u)); list(i = i, u = u, F = Fi, V = V, G = t(Fi) %*% solve(V, Fi), b = t(Fi) %*% solve(V, u)) })
G <- Reduce(`+`, lapply(pieces, `[[`, "G")); b <- Reduce(`+`, lapply(pieces, `[[`, "b"))
dmv <- function(r, S) { -.5 * (length(r) * log(2 * pi) + as.numeric(determinant(S)$modulus) + sum(r * solve(S, r))) }
chk <- sapply(pieces, function(p) { Gm <- G - p$G; bm <- b - p$b; dmv(p$u - as.vector(p$F %*% solve(Gm, bm)), p$V + p$F %*% solve(Gm, t(p$F))) })
cat("treated: max abs difference, sampler vs direct:", max(abs(chk - ll[which(trt), s])), "\n")
ctl <- which(!trt)[1:5]
chk2 <- sapply(ctl, function(i) { di <- d[d$id == i, ]; P <- pc[, round(di$t, 3) * 1000 + 1, drop = FALSE]; g <- crossprod(P, Mu) + cbind(di$X1, di$X2) %*% Bc + di$Z1 * crossprod(P, Be[, 1]) + di$Z2 * crossprod(P, Be[, 2]); Fi <- t(P) %*% Th; dmv(di$y0 - as.vector(g), Fi %*% Dm %*% t(Fi) + sg * diag(nrow(di))) })
cat("controls: max abs difference:", max(abs(chk2 - ll[ctl, s])), "\n")
x <- elpd_cut(ll, !trt); cat("elpd", round(x$elpd, 2), "ctl", round(x$elpd_ctl, 2), "trt", round(x$elpd_trt, 2), " max pareto k", round(max(x$pareto_k), 2), "\n")
# the chain itself is unchanged by the leave-one-out mode
set.seed(11); f0 <- fit_gsc(d, 2, nloop = 300, burnin = 200, psloo = 0)
cat("chain unchanged:", isTRUE(all.equal(f[[10]][[50]], f0[[10]][[50]])), isTRUE(all.equal(f[[6]][[99]], f0[[6]][[99]])), "\n")
