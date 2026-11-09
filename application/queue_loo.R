# Rerun the fits behind the k comparison with the leave-one-out calculation for the cut posterior. Same seeds as before,
# so the chains and ATT estimates should reproduce exactly; this is checked at the end of each run.
source("app2_core.R")
dir.create("res2_oldloo", showWarnings = FALSE)
long <- list(nloop = 5500, burnin = 500); short <- list(nloop = 3500, burnin = 500)
specs <- list(
  c(list(label = "k3_s1", Q = 3, seed = 20264104, psloo = 1), long), c(list(label = "k2_s1", Q = 2, seed = 20263104, psloo = 1), long),
  c(list(label = "primary_s1", Q = 1, seed = 20262104, psloo = 1), long), c(list(label = "kk4", Q = 4, seed = 20261104, psloo = 1), short),
  c(list(label = "kk5", Q = 5, seed = 20261104, psloo = 1), short), c(list(label = "kk3", Q = 3, seed = 20261104, psloo = 1), short),
  c(list(label = "k2_s2", Q = 2, seed = 20265104, psloo = 1), long))
for (sp in specs) {
  f <- file.path("res2", paste0(sp$label, ".rds")); old <- if (file.exists(f)) readRDS(f) else NULL
  if (!is.null(old)) file.copy(f, file.path("res2_oldloo", basename(f)), overwrite = TRUE)
  t0 <- Sys.time(); r <- run_app(sp)
  same <- if (is.null(old)) NA else isTRUE(all.equal(r$att$est, old$att$est, tolerance = 1e-10))
  cat(format(Sys.time(), "%H:%M:%S"), sp$label, "k", sp$Q, "secs", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "ATT unchanged", same,
      "elpd", round(r$loo$elpd, 2), "ctl", round(r$loo$elpd_ctl, 2), "trt", round(r$loo$elpd_trt, 2), "old elpd", if (!is.null(old$loo)) round(old$loo$estimates["elpd_loo", 1], 2) else NA, "max k", round(max(r$loo$pareto_k), 2), "\n")
}
L <- list(k1 = readRDS("res2/primary_s1.rds")$loo, k2 = readRDS("res2/k2_s1.rds")$loo, k3 = readRDS("res2/k3_s1.rds")$loo, k4 = readRDS("res2/kk4.rds")$loo, k5 = readRDS("res2/kk5.rds")$loo)
for (a in c("k1", "k2", "k4", "k5")) { x <- elpd_diff(L$k3, L[[a]]); ctl <- L$k3$is_control; dc <- L$k3$elpd_i[ctl] - L[[a]]$elpd_i[ctl]; cat("k3 minus", a, "diff", round(x[1], 1), "se", round(x[2], 1), " controls only", round(sum(dc), 1), "se", round(sqrt(length(dc)) * sd(dc), 1), "\n") }
cat("queue_loo finished\n")
