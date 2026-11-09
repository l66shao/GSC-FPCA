# Exact refits for controls whose Pareto k exceeds 0.7, for the fits behind the k comparison, plus the negative control with site.
source("app2_core.R")
dir.create("reloo", showWarnings = FALSE)
lme <- function(v) { m <- max(v); m + log(mean(exp(v - m))) }
fits <- list(list(label = "primary_s1", Q = 1), list(label = "k2_s1", Q = 2), list(label = "k3_s1", Q = 3), list(label = "kk4", Q = 4), list(label = "kk5", Q = 5))
sites <- c("siteB", "siteC", "siteD", "siteE")
if (!file.exists("res2/k2_negctl_icv_site.rds")) {
  r <- run_app(list(label = "k2_negctl_icv_site", Q = 2, outcome = "sri24_icv_brain_volume", covs_t = c("sexF", "race_labelBlack", "race_labelAsian", "race_labelWhite"), covs_c = sites))
  a <- r$att; print(a[a$method == "gsc1" & a$bin %in% c("all", "year 1", "year 2", "year 3"), c("bin", "n", "est", "lo", "hi")], digits = 3)
}
for (f in fits) {
  l <- readRDS(file.path("res2", paste0(f$label, ".rds")))$loo
  rows <- which(l$is_control)[l$pareto_k > .7]
  for (i in rows) {
    out <- file.path("reloo", sprintf("%s_%d.rds", f$label, i)); if (file.exists(out)) next
    t0 <- Sys.time()
    x <- run_app(list(label = paste0("reloo_", f$label), Q = f$Q, seed = 7000 + i, nloop = 3500, burnin = 500, holdout = i))
    v <- lme(x$ll); saveRDS(list(unit = i, elpd = v, psis = l$elpd_i[i]), out)
    cat(format(Sys.time(), "%H:%M:%S"), f$label, "k", f$Q, "unit", i, "pareto", round(l$pareto_k[match(i, which(l$is_control))], 2), "psis", round(l$elpd_i[i], 2), "exact", round(v, 2), "secs", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
  }
}
cat("queue_reloo finished\n")
