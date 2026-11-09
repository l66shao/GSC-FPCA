# Choice of k under the rule in Section 2.3 (simplest k within z SEs of the best), on the same datasets and chains as the
# k-choice scenarios. Rscript rule_check.R <worker 1 or 2>
library(Rcpp); library(MASS); library(splines); pdf(NULL); Sys.setenv(COMEGA = "1")
wk <- as.integer(commandArgs(TRUE)[1])
Rcpp::sourceCpp("../method/sampler_tm.cpp", cacheDir = "cpp_cache")
source("../method/mfpca.R"); source("base_functions.R"); source("sim2_functions.R"); source("../method/loo_cut.R")
params <- make_params(); cells <- read.csv("cells2.csv")
out <- sprintf("rule_check_w%d.csv", wk)
for (it in 1:100) { if (it %% 2 != wk %% 2) next
  for (cl in c(29, 30)) {
    a <- cells[cl, ]; base <- 1e6 * cl + it; set.seed(base)
    tc <- truth_curves(a$dgp, params); d <- simulate2(a$Ntr, a$Nco, a$lambda, a$n_min, a$dgp, tc, n_cap = a$n_cap)
    ctl <- !sapply(sort(unique(d$id)), function(i) d$trt[d$id == i][1]); L <- list()
    for (Q in 1:3) { set.seed(base + 10 + Q); f <- tryCatch(fit_gsc(d, Q, psloo = 1), error = function(e) NULL); if (!is.null(f)) L[[Q]] <- elpd_cut(f[[12]], ctl); rm(f) }
    e <- sapply(L, function(x) x$elpd); best <- which.max(e)
    pick <- function(z) { for (Q in 1:best) { x <- elpd_diff(L[[best]], L[[Q]]); if (Q == best || x[["diff"]] <= z * x[["se"]]) return(Q) } }
    write.table(data.frame(cell = cl, iter = it, best = best, z1 = pick(1), z2 = pick(2)), out, sep = ",", row.names = FALSE, col.names = !file.exists(out), append = TRUE)
  }
}
cat("rule check finished\n")
