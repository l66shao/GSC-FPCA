# Application analyses with k = 2 as the primary specification. Run with: Rscript queue_k2.R <worker>. Workers share the list through claim folders.
source("app2_core.R")
dir.create("claims_app_k2", showWarnings = FALSE)
sites <- c("siteB", "siteC", "siteD", "siteE")
base_t <- c("icv_baseline", "sexF", "race_labelBlack", "race_labelAsian", "race_labelWhite")
long <- list(nloop = 5500, burnin = 500)
specs <- list()
add <- function(...) specs[[length(specs) + 1]] <<- list(...)
for (s in 1:4) do.call(add, c(list(label = sprintf("k2_s%d", s), Q = 2, seed = 20261104 + 2000 * s, psloo = as.numeric(s == 1)), long))
add(label = "k2_placebo", Q = 2, placebo = TRUE)
add(label = "k2_left", Q = 2, outcome = "sf_left"); add(label = "k2_right", Q = 2, outcome = "sf_right"); add(label = "k2_orig", Q = 2, outcome = "sf_orig")
add(label = "k2_disc_itt", Q = 2, discont = "itt"); add(label = "k2_disc_censor", Q = 2, discont = "censor")
add(label = "k2_ctl_nolater", Q = 2, controls = "no_later"); add(label = "k2_ctl_zero", Q = 2, controls = "zero")
add(label = "k2_site", Q = 2, covs_c = sites)
add(label = "k2_rich", Q = 2, covs_t = c(base_t, "ses_yoe_b", "fh_alc_b", "pds_b"), covs_c = sites)
add(label = "k2_rich_sub", Q = 2, covs_t = c(base_t, "ses_yoe_b", "fh_alc_b", "pds_b", "mj_life_b", "tob_life_b"), covs_c = sites)
add(label = "k2_couse", Q = 2, covs_c = c(sites, "mj_life_t", "tob_life_t"))
add(label = "k2_win25", Q = 2, amax = 25); add(label = "k2_win25_censor", Q = 2, amax = 25, discont = "censor")
add(label = "k2_age_bins", Q = 2, age_bins = TRUE)
add(label = "k2_negctl_icv", Q = 2, outcome = "sri24_icv_brain_volume", covs_t = c("sexF", "race_labelBlack", "race_labelAsian", "race_labelWhite"))
add(label = "k2_total_gm", Q = 2, outcome = "sri24_icv_brain_gm")
add(label = "negctl_icv", outcome = "sri24_icv_brain_volume", covs_t = c("sexF", "race_labelBlack", "race_labelAsian", "race_labelWhite"))
add(label = "placebo", placebo = TRUE)  # k = 1 placebo, reported in Section 4
for (r in readRDS("rois.rds")) add(label = paste0("k2_roi_", r), Q = 2, outcome = r, nloop = 1500, burnin = 500)
for (r in 1:50) add(label = sprintf("k2_fake_%02d", r), Q = 2, fake = TRUE, seed = 9000 + r, nloop = 1000, burnin = 500)
for (sp in specs) {
  f <- file.path("res2", paste0(sp$label, ".rds"))
  if (file.exists(f)) next
  if (!dir.create(file.path("claims_app_k2", sp$label), showWarnings = FALSE)) next
  t0 <- Sys.time()
  r <- tryCatch(run_app(sp), error = function(e) { cat("ERROR", sp$label, conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) {
    a <- r$att[r$att$method == "gsc1" & r$att$bin %in% c("all", "placebo"), ]
    cat(format(Sys.time(), "%H:%M:%S"), sp$label, "N", r$meta$Ntr, r$meta$Nco, "post", r$meta$n_post, "secs", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "ATT", round(a$est, 4), "(", round(a$lo, 3), round(a$hi, 3), ")\n")
  }
}
cat("queue finished\n")
