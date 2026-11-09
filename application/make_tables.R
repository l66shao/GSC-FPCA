# Sensitivity table, ROI table, and a CSV of numbers for Section 4, from the application results.
res <- function(l) { f <- file.path("res2", paste0(l, ".rds")); if (file.exists(f)) readRDS(f) else NULL }
neg <- function(x, d = 3) { v <- sprintf(paste0("%.", d, "f"), abs(x)); ifelse(x < 0 & as.numeric(v) != 0, paste0("$-$", v), v) }
ci <- function(a, b, d = 3) sprintf("[%s,\\;%s]", neg(a, d), neg(b, d))
pick <- function(r, meth = "gsc1", bin = "all") { a <- r$att; a[a$method == meth & a$bin == bin, ] }
cell <- function(r, meth = "gsc1", bin = "all") { if (is.null(r)) return("--"); x <- pick(r, meth, bin); sprintf("%s %s", neg(x$est), ci(x$lo, x$hi)) }
cell2 <- function(r, meth = "gsc1", bin = "all") cell(r, meth, bin)
specs <- list(
  list("Primary specification ($k=2$)", "k2_s1", "primary_s1"),
  list("Prior mean zero", "k2_s1", "primary_s1", "gsc0"),
  list("\\emph{Number of FPCs}", NULL),
  list("$k=1$", "primary_s1", NULL), list("$k=3$", "k3_s1", NULL), list("$k=4$", "kk4", NULL), list("$k=5$", "kk5", NULL),
  list("\\emph{Outcome}", NULL),
  list("Left hemisphere only", "k2_left", "left"), list("Right hemisphere only", "k2_right", "right"),
  list("\\emph{Exposure definition}", NULL),
  list("Discontinuers as exposed", "k2_disc_itt", "disc_itt"), list("Discontinuers censored", "k2_disc_censor", "disc_censor"),
  list("\\emph{Controls}", NULL),
  list("Excluding later exposed", "k2_ctl_nolater", "ctl_nolater"), list("No binge drinking", "k2_ctl_zero", "ctl_zero"),
  list("\\emph{Covariates added}", NULL),
  list("Site", "k2_site", "site"), list("Site, education, family history, puberty", "k2_rich", "rich"), list("Above, baseline cannabis, tobacco", "k2_rich_sub", "rich_sub"), list("Site, cumulative cannabis, tobacco", "k2_couse", "couse"),
  list("\\emph{Age range}", NULL),
  list("Post-exposure visits to age 20.5", "k2_age_bins", NULL, "gsc1", "all age<=20.5"),
  list("Ages to 25, discontinuers dropped", "k2_win25", "win25"), list("Ages to 25, discontinuers censored", "k2_win25_censor", "win25_censor"))
rows <- character()
for (sp in specs) {
  if (is.null(sp[[2]])) { rows <- c(rows, sprintf("\\multicolumn{4}{l}{%s} \\\\", sp[[1]])); next }
  r2 <- res(sp[[2]]); r1 <- if (is.null(sp[[3]])) NULL else res(sp[[3]]); m <- if (length(sp) >= 4) sp[[4]] else "gsc1"; b <- if (length(sp) >= 5) sp[[5]] else "all"
  nn <- if (is.null(r2)) "--" else sprintf("%d/%d/%d", r2$meta$Ntr, r2$meta$Nco, pick(r2, m, b)$n)
  rows <- c(rows, sprintf("%s & %s & %s & %s \\\\", sp[[1]], nn, cell(r2, m, b), if (is.null(sp[[3]])) "--" else cell(r1, m, "all")))
}
af <- read.csv("out_final_k2/att_table.csv"); xf <- af[af$label == "all post-exposure", ]
rows[1] <- sprintf("Primary specification ($k=2$) & 115/513/182 & %s %s & %s \\\\", neg(xf$est), ci(xf$pred_lo, xf$pred_hi), cell(res("primary_s1")))
writeLines(rows, "tab_sens_rows.tex")
# ROI table
rois <- readRDS("rois.rds")
tab <- do.call(rbind, lapply(rois, function(ro) { r <- if (ro == "superiorfrontal_grayvol") res("k2_s1") else res(paste0("k2_roi_", ro)); if (is.null(r)) return(NULL); x <- pick(r); data.frame(roi = ro, est = x$est, lo = x$lo, hi = x$hi, p_neg = x$p_neg, sc = r$meta$sc) }))
if (!is.null(tab)) {
  af <- read.csv("out_final_k2/att_table.csv"); xf <- af[af$label == "all post-exposure", ]; js <- tab$roi == "superiorfrontal_grayvol"
  tab$est[js] <- xf$est; tab$lo[js] <- xf$pred_lo; tab$hi[js] <- xf$pred_hi; tab$p_neg[js] <- xf$post_prob_neg
  tab$p2 <- 2 * pmin(tab$p_neg, 1 - tab$p_neg); tab$p2[tab$p2 < 2 / 1000] <- 2 / 1000
  tab$bh <- p.adjust(tab$p2, "BH")
  tab <- tab[order(tab$est), ]; tab$rank <- seq_len(nrow(tab))
  write.csv(tab, "roi_table.csv", row.names = FALSE)
  nm <- function(x) { x <- sub("_grayvol$", "", x); x }
  pretty <- c(bankssts = "Banks of superior temporal sulcus", caudalanteriorcingulate = "Caudal anterior cingulate", caudalmiddlefrontal = "Caudal middle frontal", cuneus = "Cuneus", entorhinal = "Entorhinal", fusiform = "Fusiform", inferiorparietal = "Inferior parietal", inferiortemporal = "Inferior temporal", isthmuscingulate = "Isthmus cingulate", lateraloccipital = "Lateral occipital", lateralorbitofrontal = "Lateral orbitofrontal", lingual = "Lingual", medialorbitofrontal = "Medial orbitofrontal", middletemporal = "Middle temporal", parahippocampal = "Parahippocampal", paracentral = "Paracentral", parsopercularis = "Pars opercularis", parsorbitalis = "Pars orbitalis", parstriangularis = "Pars triangularis", pericalcarine = "Pericalcarine", postcentral = "Postcentral", posteriorcingulate = "Posterior cingulate", precentral = "Precentral", precuneus = "Precuneus", rostralanteriorcingulate = "Rostral anterior cingulate", rostralmiddlefrontal = "Rostral middle frontal", superiorfrontal = "Superior frontal", superiorparietal = "Superior parietal", superiortemporal = "Superior temporal", supramarginal = "Supramarginal", frontalpole = "Frontal pole", temporalpole = "Temporal pole", transversetemporal = "Transverse temporal", insula = "Insula")
  writeLines(sprintf("%s & %s & %s & %s \\\\", pretty[nm(tab$roi)], neg(tab$est), ci(tab$lo, tab$hi), sprintf("%.3f", tab$bh)), "tab_roi_rows.tex")
}
# calibration with pseudo-exposures
fk <- list.files("res2", pattern = "^k2_fake_.*rds$", full.names = TRUE)
if (length(fk)) {
  z <- do.call(rbind, lapply(fk, function(f) { x <- pick(readRDS(f)); data.frame(est = x$est, lo = x$lo, hi = x$hi, n = x$n, units = x$units) }))
  z0 <- do.call(rbind, lapply(fk, function(f) { x <- pick(readRDS(f), "gsc0"); data.frame(est = x$est) }))
  cat(sprintf("fake runs %d mean %.4f sd %.4f min %.3f max %.3f cover %d avgwidth %.3f mean_n %.1f units %.1f | gsc0 mean %.4f\n", nrow(z), mean(z$est), sd(z$est), min(z$est), max(z$est), sum(z$lo <= 0 & z$hi >= 0), mean(z$hi - z$lo), mean(z$n), mean(z$units), mean(z0$est)))
  write.csv(z, "fake_table.csv", row.names = FALSE)
}
