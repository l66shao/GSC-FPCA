library(coda)
ch <- lapply(1:4, function(s) readRDS(sprintf("res2/k2_s%d.rds", s)))
get <- function(fun) mcmc.list(lapply(ch, function(r) mcmc(fun(r))))
items <- list("ATT all" = function(r) r$draws[["gsc1 all"]], "ATT year 1" = function(r) r$draws[["gsc1 year 1"]], "ATT year 2" = function(r) r$draws[["gsc1 year 2"]], "ATT year 3" = function(r) r$draws[["gsc1 year 3"]],
 "ATT all (prior mean 0)" = function(r) r$draws[["gsc0 all"]], "sigma2" = function(r) r$sigma2)
for (j in 1:7) local({ jj <- j; items[[paste0("f'nu at grid ", jj)]] <<- function(r) r$NuF[, jj] })
for (j in 1:10) local({ jj <- j; items[[paste0("cov ", jj)]] <<- function(r) r$Cov[, jj] })
out <- do.call(rbind, lapply(names(items), function(n) { m <- get(items[[n]]); g <- gelman.diag(m, autoburnin = FALSE, multivariate = FALSE)$psrf
  data.frame(item = n, mean = mean(unlist(m)), sd = sd(unlist(m)), rhat = g[1], rhat_u = g[2], ess_total = sum(effectiveSize(m)), ess_min_chain = min(sapply(m, effectiveSize))) }))
print(out, digits = 4, row.names = FALSE)
write.csv(out, "diag_k2.csv", row.names = FALSE)
est <- sapply(ch, function(r) r$att$est[r$att$method == "gsc1" & r$att$bin == "all"]); cat("chain ATT all:", round(est, 4), " range", round(diff(range(est)), 4), "\n")
cat("age grid", round(ch[[1]]$age_grid, 1), "\n"); nf <- do.call(rbind, lapply(ch, function(r) r$NuF)); print(round(rbind(mean = colMeans(nf), lo = apply(nf, 2, quantile, .025), hi = apply(nf, 2, quantile, .975)), 3))
