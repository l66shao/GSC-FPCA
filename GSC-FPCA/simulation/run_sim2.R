# GSC-FPCA expanded simulation. Each worker claims tasks (cell, replicate) from a shared queue, so workers can be added,
# removed or restarted at any time: Rscript run_sim2.R <worker>. Results are appended by task; a task is finished when
# claims/<cell>_<iter>.done exists. After a crash, with no workers running, remove claim folders lacking a .done file.
#
# Scenarios. Exposure and sparsity are as in Section 3 unless noted; the effect is 5*delta, so errors of estimators that
# impute from untreated cells do not depend on its shape. Outcomes in "none" follow the Section 3 model.
#   none    exposure unrelated to scores
#   select  the Ntr largest of U = alpha_1 + alpha_2 + logistic noise are exposed
#   growth  intercept and slope (correlation .3) rather than FPCs, Ornstein-Uhlenbeck plus white errors, selection as above
#   outside mean and FPCs outside the fitted spline span, Gamma(2) scores, t_5 errors, selection as above
#   wrongk  three components, two fitted, selection as above
#   visits  visit counts depend on the first score, selection as above
#   dip     exposure comes right after a visit with an unusually low error, selection as above; no method should cope
#   dense   every unit observed at all 20 grid points, selection as above
#   looic   k chosen from {1,2,3} by LOOIC under select (true k=2) and wrongk (true k=3)
# Methods: gsc1 (separate prior mean for exposed units), gsc0 (prior mean zero), did (imputation), lme1, lme2 (mixed models),
#          ife2 (interactive fixed effects with r=2), ifecv (r by cross validation), mc (matrix completion).
library(Rcpp); library(MASS); library(splines); pdf(NULL); Sys.setenv(COMEGA = "1")
worker <- commandArgs(TRUE)[1]
dir.create("results2", showWarnings = FALSE); dir.create("claims", showWarnings = FALSE); dir.create("cpp_cache", showWarnings = FALSE)
Rcpp::sourceCpp("../method/sampler_tm.cpp", cacheDir = "cpp_cache")
source("../method/mfpca.R"); source("base_functions.R"); source("sim2_functions.R"); source("../method/loo_cut.R")
params <- make_params(); Bm <- truth_basis()
bins <- c(sprintf("%.1f", (1:9) / 10), "all")
sp <- rbind(A = c(2, 2), B = c(4, 3), C = c(7, 4)); sz <- rbind(s = c(20, 50), m = c(40, 100), l = c(80, 200))
cells <- list(); add <- function(scen, dgp, s, p, reps, r_lme = 0, r_ife2 = 0, r_cv = 0, r_mc = 0, r_shape = 0, r_boot = 0, Q = 2, n_cap = 10, lambda = NULL, n_min = NULL, looic = FALSE) {
  cells[[length(cells) + 1]] <<- data.frame(scen = scen, dgp = dgp, size = s, sp = p, Ntr = sz[s, 1], Nco = sz[s, 2], lambda = if (is.null(lambda)) sp[p, 1] else lambda, n_min = if (is.null(n_min)) sp[p, 2] else n_min,
    n_cap = n_cap, reps = reps, r_lme = r_lme, r_ife2 = r_ife2, r_cv = r_cv, r_mc = r_mc, r_shape = r_shape, r_boot = r_boot, Q = Q, looic = looic)
}
for (scen in c("none", "select")) for (s in c("s", "m", "l")) for (p in c("A", "B", "C")) {
  if (s == "m") add(scen, scen, s, p, 300, r_lme = 200, r_ife2 = 100, r_cv = 100, r_mc = 200, r_shape = if (p == "B") 100 else 0, r_boot = if (p == "B") 60 else 0)
  else add(scen, scen, s, p, 300)
}
for (scen in c("growth", "outside")) for (p in c("A", "B", "C")) add(scen, scen, "m", p, 200, r_lme = 200, r_ife2 = 100, r_cv = 100, r_mc = 200)
for (scen in c("wrongk", "visits", "dip")) add(scen, scen, "m", "B", 200, r_lme = 200, r_ife2 = 100, r_cv = 100, r_mc = 200)
add("dense", "select", "m", "C", 100, r_lme = 100, r_ife2 = 100, r_cv = 100, r_mc = 100, n_cap = 20, lambda = 100, n_min = 20)
add("looic", "select", "m", "B", 200, looic = TRUE); add("looic3", "wrongk", "m", "B", 200, looic = TRUE)
cells <- do.call(rbind, cells); cells$cell <- seq_len(nrow(cells)); rownames(cells) <- NULL
if (!file.exists("cells2.csv")) write.csv(cells, "cells2.csv", row.names = FALSE)
tasks <- do.call(rbind, lapply(1:max(cells$reps), function(it) { w <- which(cells$reps >= it); data.frame(cell = w, iter = it) }))
if (!is.na(Sys.getenv("SIM_MAX", NA))) tasks <- tasks[tasks$iter <= as.numeric(Sys.getenv("SIM_MAX")), ]
if (!is.na(Sys.getenv("SIM_CELLS", NA))) tasks <- tasks[tasks$cell %in% eval(parse(text = Sys.getenv("SIM_CELLS"))), ]
outfile <- sprintf("results2/w%s.csv", worker)
# release this worker's unfinished claims from an earlier run
for (cdir in list.files("claims", full.names = TRUE)) {
  if (dir.exists(cdir) && file.exists(file.path(cdir, paste0("w", worker))) && !file.exists(paste0(cdir, ".done"))) unlink(cdir, recursive = TRUE)
}
cols <- c("cell", "iter", "scen", "size", "sp", "method", "shape", "bin", "n", "err", "lo", "hi", "indiv_mse", "sigma2", "kept", "r_used", "truth_linear", "chosen_k", "nu1", "nu2", "secs")
if (!file.exists(outfile)) write.table(t(cols), outfile, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
safely <- function(expr, lab) tryCatch(expr, error = function(e) { cat(sprintf("ERROR %s: %s\n", lab, conditionMessage(e))); NULL })
timed <- function(expr) { t0 <- Sys.time(); r <- force(expr); list(r = r, secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))) }
for (ti in seq_len(nrow(tasks))) {
  cl <- tasks$cell[ti]; it <- tasks$iter[ti]
  claim <- sprintf("claims/%d_%d", cl, it)
  if (!dir.create(claim, showWarnings = FALSE)) next
  file.create(file.path(claim, paste0("w", worker)))
  a <- cells[cl, ]
  base <- 1e6 * cl + it
  set.seed(base)
  tc <- truth_curves(a$dgp, params)
  d <- simulate2(a$Ntr, a$Nco, a$lambda, a$n_min, a$dgp, tc, n_cap = a$n_cap)
  d$bin <- ifelse(d$post, sprintf("%.1f", round(d$delta * 10) / 10), NA)
  nu <- attr(d, "nu")
  res <- list()
  keep <- function(x, secs, lab = "") {
    if (is.null(x)) { res[[length(res) + 1]] <<- data.frame(method = lab, shape = "any", bin = "all", n = NA, err = NA, lo = NA, hi = NA, secs = secs); return(invisible()) }
    if (!"shape" %in% names(x)) x$shape <- "any"
    x$secs <- secs
    res[[length(res) + 1]] <<- x
  }
  if (a$looic) {
    elpd <- rep(NA, 3); gs <- list()
    for (Q in 1:3) {
      set.seed(base + 10 + Q)
      tf <- timed(safely(fit_gsc(d, Q, psloo = 1), sprintf("gsc k=%d cell %d iter %d", Q, cl, it)))
      if (is.null(tf$r)) next
      elpd[Q] <- elpd_cut(tf$r[[12]], !sapply(sort(unique(d$id)), function(i) d$trt[d$id == i][1]))$elpd
      gs[[Q]] <- gsc_att(tf$r, d, bins); gs[[Q]]$secs <- tf$secs; gs[[Q]]$r_used <- Q
      rm(tf)
    }
    ch <- which.max(elpd)
    for (Q in which(!is.na(elpd))) { x <- gs[[Q]]; x$chosen_k <- ch; x$elpd <- elpd[Q]; x$method <- paste0(x$method, "_k", Q); x$elpd <- NULL; keep(x, x$secs[1]) }
    x <- gs[[ch]]; x$chosen_k <- ch; x$method <- paste0(x$method, "_looic"); keep(x, 0)
  } else {
    set.seed(base + 1)
    tf <- timed(safely(fit_gsc(d, a$Q), sprintf("gsc cell %d iter %d", cl, it)))
    if (!is.null(tf$r)) { x <- gsc_att(tf$r, d, bins); keep(x, tf$secs / 2, "gsc") } else keep(NULL, tf$secs, "gsc1")
    rm(tf)
    set.seed(base + 2)
    tf <- timed(safely(did_att(d, Bm, bins), "did")); keep(tf$r, tf$secs, "did")
    if (it <= a$r_lme) {
      set.seed(base + 3)
      shp <- if (it <= a$r_shape) c("linear", "constant", "saturating", "transient") else "linear"
      tf <- timed(safely(lme_att(d, Bm, bins, shp), "lme")); keep(tf$r, tf$secs, "lme1")
    }
    if (it <= a$r_ife2) { set.seed(base + 4); tf <- timed(safely(fect_att(d, Bm, bins, method = "ife", r = 2, min_T0 = 3, label = "ife2"), "ife2")); keep(tf$r, tf$secs, "ife2") }
    if (it <= a$r_cv) { set.seed(base + 5); tf <- timed(safely(fect_att(d, Bm, bins, method = "ife", cv = TRUE, min_T0 = 3, label = "ifecv"), "ifecv")); keep(tf$r, tf$secs, "ifecv") }
    if (it <= a$r_mc) { set.seed(base + 6); tf <- timed(safely(fect_att(d, Bm, bins, method = "mc", min_T0 = 1, se = it <= a$r_boot, nboots = 100, label = "mc"), "mc")); keep(tf$r, tf$secs, "mc") }
  }
  r <- do.call(rbind, lapply(res, function(x) {
    for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
    x$cell <- cl; x$iter <- it; x$scen <- a$scen; x$size <- a$size; x$sp <- a$sp; x$nu1 <- nu[1]; x$nu2 <- nu[2]
    x[, cols]
  }))
  # round to keep the files small
  for (nm in c("err", "lo", "hi", "indiv_mse", "sigma2", "truth_linear", "nu1", "nu2", "secs", "kept")) r[[nm]] <- round(r[[nm]], 5)
  write.table(r, outfile, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
  file.create(paste0(claim, ".done"))
  cat(format(Sys.time(), "%H:%M:%S"), "cell", cl, a$scen, a$size, a$sp, "iter", it, "secs", round(sum(sapply(res, function(x) x$secs[1]), na.rm = TRUE), 1), "\n")
}
cat("worker", worker, "finished\n")
