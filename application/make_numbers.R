# Numbers for Section 4 and Appendix B, written to tokens_app.json. Primary numbers come from the chain that produced the figures.
library(jsonlite)
OUT <- "out_final_k2"; RES <- "res2"
res <- function(l) readRDS(file.path(RES, paste0(l, ".rds")))
neg <- function(x, d = 3) { v <- sprintf(paste0("%.", d, "f"), abs(x)); ifelse(x < 0 & as.numeric(v) != 0, paste0("$-$", v), v) }
ci <- function(a, b, d = 3) sprintf("$[%s,\\;%s]$", gsub("\\$", "", neg(a, d)), gsub("\\$", "", neg(b, d)))
mneg <- function(x, d = 3) sprintf("$%s$", gsub("\\$", "", neg(x, d)))
pick <- function(r, meth = "gsc1", bin = "all") { a <- r$att; a[a$method == meth & a$bin == bin, ] }
a <- read.csv(file.path(OUT, "att_table.csv")); sc <- a$sc[1]
T <- list()
row <- function(l) a[a$label == l, ]
for (k in 1:3) { x <- row(sprintf("year %d", k)); T[[paste0("n", k)]] <- x$n; T[[paste0("e", k)]] <- mneg(x$est); T[[paste0("i", k)]] <- ci(x$pred_lo, x$pred_hi); T[[paste0("o", k)]] <- mneg(x$est * sc, 1); T[[paste0("oi", k)]] <- ci(x$pred_lo * sc, x$pred_hi * sc, 1) }
x <- row("all post-exposure"); T$nall <- x$n; T$eall <- mneg(x$est); T$iall <- ci(x$pred_lo, x$pred_hi); T$oall <- mneg(x$est * sc, 1); T$oiall <- ci(x$pred_lo * sc, x$pred_hi * sc, 1)
# ATT table in original units for Appendix B
rw <- character()
for (k in -3:0) { x <- row(sprintf("year %d", k)); rw <- c(rw, sprintf("$%s$ & $%s$ & $[%s,\\;%s]$ & $%d$ \\\\", gsub("\\$", "", neg(k, 0)), gsub("\\$", "", neg(x$est * sc, 1)), gsub("\\$", "", neg(x$fit_lo * sc, 1)), gsub("\\$", "", neg(x$fit_hi * sc, 1)), x$n)) }
rw <- c(rw, "\\hline")
for (k in 1:3) { x <- row(sprintf("year %d", k)); rw <- c(rw, sprintf("$%d$ & $%s$ & $[%s,\\;%s]$ & $%d$ \\\\", k, gsub("\\$", "", neg(x$est * sc, 1)), gsub("\\$", "", neg(x$pred_lo * sc, 1)), gsub("\\$", "", neg(x$pred_hi * sc, 1)), x$n)) }
writeLines(rw, "tab_att_rows.tex")
# LOO
source("../method/loo_cut.R")
l <- readRDS(file.path(OUT, "loo.rds")); pk <- l$pareto_k; N <- length(pk); T$k_n <- N
T$elpd <- gsub("\\$", "", neg(l$elpd, 1)); T$elpd_se <- sprintf("%.1f", l$se); T$looic <- sprintf("%.1f", l$looic); T$looic_se <- sprintf("%.1f", l$looic_se); T$ploo <- sprintf("%.1f", l$p_loo_ctl); T$ploo_se <- sprintf("%.1f", l$p_loo_ctl_se)
ok <- sum(pk <= .7); mid <- sum(pk > .7 & pk <= 1); hi <- sum(pk > 1)
T$k_ok <- ok; T$k_ok_pct <- sprintf("%.1f", 100 * ok / N); T$k_mid <- mid; T$k_mid_pct <- sprintf("%.1f", 100 * mid / N); T$k_hi <- hi; T$k_hi_pct <- sprintf("%.1f", 100 * hi / N)
T$k_over_word <- c("No", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten")[mid + hi + 1]
# Replace PSIS by the exact refit for controls with a refit available for that k
exact <- function(x, label) { fs <- list.files("reloo", pattern = paste0("^", label, "_[0-9]+[.]rds$"), full.names = TRUE)
  for (f in fs) { z <- readRDS(f); x$p_loo_ctl <- x$p_loo_ctl + (x$elpd_i[z$unit] - z$elpd); x$elpd_i[z$unit] <- z$elpd }
  n <- length(x$elpd_i); x$elpd <- sum(x$elpd_i); x$se <- sqrt(n) * sd(x$elpd_i); x$looic <- -2 * x$elpd; x$looic_se <- 2 * x$se; x$elpd_ctl <- sum(x$elpd_i[x$is_control]); x$elpd_trt <- sum(x$elpd_i[!x$is_control]); x$n_exact <- length(fs); x }
l <- exact(l, "k2_s1")
T$elpd <- gsub("\\$", "", neg(l$elpd, 1)); T$elpd_se <- sprintf("%.1f", l$se); T$looic <- sprintf("%.1f", l$looic); T$looic_se <- sprintf("%.1f", l$looic_se); T$ploo <- sprintf("%.1f", l$p_loo_ctl)
L <- list(k1 = exact(res("primary_s1")$loo, "primary_s1"), k2 = l, k3 = exact(res("k3_s1")$loo, "k3_s1"), k4 = exact(res("kk4")$loo, "kk4"), k5 = exact(res("kk5")$loo, "kk5"))
cat("refits used:", sapply(L, function(x) x$n_exact), "\n")
K2 <- list(l, exact(res("k2_s1")$loo, "k2_s1"), exact(res("k2_s2")$loo, "k2_s1")); K3 <- list(exact(res("k3_s1")$loo, "k3_s1"), exact(res("kk3")$loo, "k3_s1"))
e2 <- sapply(K2, function(x) x$elpd); e3 <- sapply(K3, function(x) x$elpd); cat("elpd k=2 across chains:", round(e2, 2), " k=3:", round(e3, 2), "\n")
T$elpd2_lo <- mneg(min(e2), 1); T$elpd2_hi <- mneg(max(e2), 1)
print(round(sapply(L, function(x) c(elpd = x$elpd, se = x$se, ctl = x$elpd_ctl, trt = x$elpd_trt, ploo_ctl = x$p_loo_ctl)), 2))
# k = 3 against k = 2 over all pairs of chains, with and without controls whose Pareto k exceeds 0.7 in either fit
pair <- function(a, b) { d <- a$elpd_i - b$elpd_i; ct <- a$is_control; bad <- rep(FALSE, length(d)); bad[ct] <- a$pareto_k > .7 | b$pareto_k > .7; c(d = sum(d), se = sqrt(length(d)) * sd(d), dx = sum(d[!bad]), sex = sqrt(sum(!bad)) * sd(d[!bad]), nbad = sum(bad), dt = sum(d[!ct]), set = sqrt(sum(!ct)) * sd(d[!ct])) }
P32 <- do.call(rbind, lapply(K3, function(a) do.call(rbind, lapply(K2, function(b) pair(a, b))))); print(round(P32, 2))
T$e32_lo <- sprintf("%.1f", min(P32[, "d"])); T$e32_hi <- sprintf("%.1f", max(P32[, "d"])); T$e32_se <- sprintf("%.1f", mean(P32[, "se"]))
T$e32_rlo <- sprintf("%.1f", min(P32[, "d"] / P32[, "se"])); T$e32_rhi <- sprintf("%.1f", max(P32[, "d"] / P32[, "se"]))
T$e32x_lo <- sprintf("%.1f", min(P32[, "dx"])); T$e32x_hi <- sprintf("%.1f", max(P32[, "dx"])); T$e32x_se <- sprintf("%.1f", mean(P32[, "sex"])); T$nbad_lo <- min(P32[, "nbad"]); T$nbad_hi <- max(P32[, "nbad"])
P21 <- do.call(rbind, lapply(K2, function(a) pair(a, L$k1))); print(round(P21, 2))
T$e21_lo <- sprintf("%.1f", min(P21[, "d"])); T$e21_hi <- sprintf("%.1f", max(P21[, "d"])); T$e21_se <- sprintf("%.1f", mean(P21[, "se"]))
T$e21x_lo <- sprintf("%.1f", min(P21[, "dx"])); T$e21x_hi <- sprintf("%.1f", max(P21[, "dx"])); T$e21x_se <- sprintf("%.1f", mean(P21[, "sex"]))
T$e21t <- mneg(mean(P21[, "dt"]), 1); T$e21t_se <- sprintf("%.1f", mean(P21[, "set"])); T$e32t <- mneg(mean(P32[, "dt"]), 1); T$e32t_se <- sprintf("%.1f", mean(P32[, "set"]))
P31 <- do.call(rbind, lapply(K3, function(a) pair(a, L$k1))); print(round(P31, 2)); T$e31_lo <- sprintf("%.1f", min(P31[, "d"])); T$e31_hi <- sprintf("%.1f", max(P31[, "d"])); T$e31_se <- sprintf("%.1f", mean(P31[, "se"]))
P43 <- do.call(rbind, lapply(K3, function(b) pair(L$k4, b))); P53 <- do.call(rbind, lapply(K3, function(b) pair(L$k5, b))); print(round(rbind(P43, P53), 2))
T$e43_lo <- mneg(min(P43[, "d"]), 1); T$e43_hi <- mneg(max(P43[, "d"]), 1); T$e43_se <- sprintf("%.1f", mean(P43[, "se"])); T$e53_lo <- mneg(min(P53[, "d"]), 1); T$e53_hi <- mneg(max(P53[, "d"]), 1); T$e53_se <- sprintf("%.1f", mean(P53[, "se"]))
# diagnostics
dg <- read.csv("diag_k2.csv"); T$rhat <- sprintf("%.3f", max(dg$rhat)); T$ess <- format(round(min(dg$ess_total) / 100) * 100, big.mark = ",")
ests <- sapply(1:4, function(s) pick(res(sprintf("k2_s%d", s)))$est); T$chain_range <- sprintf("%.3f", diff(range(c(ests, a$est[a$label == "all post-exposure"])))); print(ests)
# placebo
p <- pick(res("k2_placebo"), "gsc1", "placebo"); T$pl2 <- mneg(p$est); T$pl2i <- ci(p$lo, p$hi); T$pl2o <- mneg(p$est * sc, 1); T$pl2oi <- ci(p$lo * sc, p$hi * sc, 1)
p1 <- pick(res("placebo"), "gsc1", "placebo"); T$pl1 <- mneg(p1$est); T$pl1i <- ci(p1$lo, p1$hi)
# other k
for (nm in list(c("k1", "primary_s1"), c("k3", "k3_s1"), c("k4", "kk4"), c("k5", "kk5"))) { x <- pick(res(nm[2])); T[[paste0(nm[1], "all")]] <- mneg(x$est); T[[paste0(nm[1], "alli")]] <- ci(x$lo, x$hi) }
r1 <- res("primary_s1"); T$k1r3 <- mneg(pick(r1, "gsc1", "pre -3")$est); T$k1r0 <- mneg(pick(r1, "gsc1", "pre 0")$est)
T$kpr3 <- mneg(row("year -3")$est); T$kpr0 <- mneg(row("year 0")$est)
ef1 <- abs(res("primary_s1")$eig_fun[, 1]); T$f1_lo <- sprintf("%.2f", min(ef1)); T$f1_hi <- sprintf("%.2f", max(ef1))
# treated mean on the outcome scale, pooled over the chains
nf <- do.call(rbind, lapply(1:4, function(s) res(sprintf("k2_s%d", s))$NuF)); ag <- res("k2_s1")$age_grid
print(round(rbind(age = ag, mean = colMeans(nf), lo = apply(nf, 2, quantile, .025), hi = apply(nf, 2, quantile, .975)), 3))
T$nuf_12 <- mneg(mean(nf[, 1])); T$nuf_21 <- mneg(mean(nf[, 7])); T$nuf_12i <- ci(quantile(nf[, 1], .025), quantile(nf[, 1], .975)); T$nuf_21i <- ci(quantile(nf[, 7], .025), quantile(nf[, 7], .975))
g0 <- pick(res("k2_s1"), "gsc0"); T$g0all <- mneg(g0$est); T$g0alli <- ci(g0$lo, g0$hi)
# negative control and ROI
ic <- pick(res("k2_negctl_icv")); T$icv <- mneg(ic$est); T$icvi <- ci(ic$lo, ic$hi); ic1 <- pick(res("negctl_icv")); T$icv_k1 <- mneg(ic1$est); T$icv_k1i <- ci(ic1$lo, ic1$hi); ic3 <- pick(res("k3_negctl_icv")); T$icv_k3 <- mneg(ic3$est); T$icv_k3i <- ci(ic3$lo, ic3$hi); tg <- pick(res("k2_total_gm")); T$tgm <- mneg(tg$est); T$tgmi <- ci(tg$lo, tg$hi)
if (file.exists("roi_table.csv")) { R <- read.csv("roi_table.csv"); T$roi_neg <- sum(R$est < 0); T$roi_pos_name <- R$roi[which.max(R$est)]; T$roi_excl <- sum(R$hi < 0 | R$lo > 0); T$roi_bh05 <- sum(R$bh < .05); T$roi_sf_rank <- R$rank[R$roi == "superiorfrontal_grayvol"]; print(R[order(R$est)[1:8], c("roi", "est", "lo", "hi", "bh")], digits = 3); print(R[R$roi == "superiorfrontal_grayvol", ], digits = 3); T$roi_median <- mneg(median(R$est)) ; T$roi_min <- mneg(min(R$est)); T$roi_max <- mneg(max(R$est)) }
if (file.exists("fake_table.csv")) { Z <- read.csv("fake_table.csv"); T$fk_mean <- mneg(mean(Z$est)); T$fk_sd <- sprintf("%.3f", sd(Z$est)); T$fk_cov <- sum(Z$lo <= 0 & Z$hi >= 0); T$fk_n <- nrow(Z); T$fk_width <- sprintf("%.3f", mean(Z$hi - Z$lo)) }
# LME
m <- read.csv(file.path(OUT, "lme.csv")); rownames(m) <- m$X
T$lme_est <- mneg(m["binge12", "Value"]); T$lme_se <- sprintf("%.3f", m["binge12", "Std.Error"]); pv <- m["binge12", "p.value"]; T$lme_p <- sprintf("%.1f\\times 10^{%d}", pv / 10^floor(log10(pv)), floor(log10(pv)))
fp <- function(p) if (p < .001) "$<0.001$" else sprintf("%.3f", p)
lab <- c("(Intercept)" = "(Intercept)", binge12 = "Heavy binge drinking (12+ episodes/year)", icv_baseline = "Baseline intracranial volume (standardized)", sexM = "Male (vs female)", visit = "Visit (year)", base_age = "Baseline age (years)", "race_labelAfrican-American/Black" = "Race: African-American/Black", race_labelAsian = "Race: Asian", "race_labelCaucasian/White" = "Race: Caucasian/White")
writeLines(sapply(names(lab), function(n) sprintf("%s & %s & %.3f & %s \\\\", lab[[n]], mneg(m[n, "Value"]), m[n, "Std.Error"], fp(m[n, "p.value"]))), "tab_lme_app_rows.tex")
# Covariate text inputs
print(read.csv(file.path(OUT, "coef_bands.csv")), digits = 3)
print(read.csv(file.path(OUT, "fpc_eigen.csv")), digits = 3)
xa <- row("all post-exposure"); T$eall_abs <- sprintf("%.3f", abs(xa$est)); T$iall_abs <- sprintf("%.3f to %.3f", abs(xa$pred_hi), abs(xa$pred_lo))
if (file.exists("roi_table.csv")) T$roi_sf_word <- c("largest", "second largest", "third largest", "fourth largest", "fifth largest", "sixth largest", "seventh largest", "eighth largest", "ninth largest", "tenth largest", "11th largest", "12th largest")[T$roi_sf_rank]
mm3 <- function(x) gsub(",", "{,}", sprintf("$%s$", format(round(x), big.mark = ",")))
ric <- res("k2_negctl_icv"); ai <- pick(ric); T$icv_mm3 <- mm3(ai$est * ric$meta$sc); y3 <- pick(ric, "gsc1", "year 3"); T$icv_y3 <- mneg(y3$est); T$icv_y3i <- ci(y3$lo, y3$hi)
rs <- res("k2_negctl_icv_site"); as <- pick(rs); T$icv_site <- mneg(as$est); T$icv_sitei <- ci(as$lo, as$hi)
rtg <- res("k2_total_gm"); T$tgm_mm3 <- mm3(pick(rtg)$est * rtg$meta$sc)
est0 <- row("all post-exposure")$est; T$disc_itt_pct <- sprintf("%.0f", 100 * (1 - pick(res("k2_disc_itt"))$est / est0)); T$disc_cens_pct <- sprintf("%.0f", 100 * (1 - pick(res("k2_disc_censor"))$est / est0))
# contrasts between years, with the measurement-error draws that the predictive intervals include
ch2 <- lapply(1:4, function(s) res(sprintf("k2_s%d", s))); set.seed(20261105)
yb <- function(b, n) unlist(lapply(ch2, function(r) { x <- r$draws[[paste("gsc1", b)]]; x + rnorm(length(x), 0, sqrt(r$sigma2 / n)) }))
Y1 <- yb("year 1", 83); Y2 <- yb("year 2", 62); Y3 <- yb("year 3", 21)
T$d32 <- mneg(mean(Y3 - Y2)); T$d32i <- ci(quantile(Y3 - Y2, .025), quantile(Y3 - Y2, .975)); T$d31 <- mneg(mean(Y3 - Y1)); T$d31i <- ci(quantile(Y3 - Y1, .025), quantile(Y3 - Y1, .975))
write_json(T, "tokens_app.json", auto_unbox = TRUE, pretty = TRUE)
