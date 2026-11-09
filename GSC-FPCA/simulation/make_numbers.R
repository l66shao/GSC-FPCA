# Numbers quoted in the prose of Section 3, read from summary2.csv so that the text follows the tables.
library(data.table); library(jsonlite)
s <- fread("summary2.csv"); ch <- fread("looic_choice.csv")
f <- function(x, d = 2) { v <- sprintf(paste0("%.", d, "f"), abs(x)); ifelse(x < 0 & as.numeric(v) != 0, paste0("$-$", v), v) }
fm <- function(x, d = 2) { v <- sprintf(paste0("%.", d, "f"), abs(x)); ifelse(x < 0 & as.numeric(v) != 0, paste0("-", v), v) }   # for use inside math mode
num <- function(x, d = 2) sprintf(paste0("%.", d, "f"), x)
main <- s[bin == "all" & scen %in% c("none", "select") & shape == "any"]
g <- function(sc, m) main[scen == sc & method == m]
T <- list()
ng <- function(sc, m = "gsc1", extra = TRUE) min(s[scen %in% sc & method == m & bin == "all" & shape %in% c("any", "linear")]$reps)
R <- ng(c("none", "select")); T$R <- as.character(R); T$mcse_cov <- num(sqrt(.95 * .05 / R), 3)
T$R_other <- as.character(ng(c("growth", "outside", "wrongk", "visits", "dip"))); T$R_dense <- as.character(ng("dense"))
loo_x <- s[scen %in% c("looic", "looic3") & method == "gsc1_looic" & bin == "all"]; T$R_loo <- as.character(min(loo_x$reps))
sm <- s[size == "m" & sp == "B" & scen %in% c("none", "select") & bin == "all"]
T$R_lme <- as.character(min(sm[method == "lme1" & shape == "linear"]$reps)); T$R_mc <- as.character(min(sm[method == "mc"]$reps)); T$R_ife <- as.character(min(sm[method == "ife2"]$reps)); T$R_cv <- as.character(min(sm[method == "ifecv"]$reps))
T$rb <- as.character(s[scen == "none" & size == "m" & sp == "B" & method == "mc" & bin == "all"]$ncov)
# scores
x <- s[scen == "select" & method == "gsc1" & bin == "all"]
T$nu1 <- num(mean(x$nu1), 2); T$nu2 <- num(mean(x$nu2), 2)
# no selection / selection, proposed estimator
for (sc in c("none", "select")) {
  x <- g(sc, "gsc1"); p <- if (sc == "none") "n" else "s"
  T[[paste0(p, "_maxb")]] <- num(max(abs(x$bias)), 3)
  T[[paste0(p, "_n2se")]] <- as.character(sum(abs(x$bias) > 2 * x$se_bias))
  T[[paste0(p, "_covlo")]] <- num(min(x$cov), 3); T[[paste0(p, "_covhi")]] <- num(max(x$cov), 3)
  T[[paste0(p, "_wlo")]] <- num(min(x$width), 2); T[[paste0(p, "_whi")]] <- num(max(x$width), 2)
  T[[paste0(p, "_sdSA")]] <- num(x[size == "s" & sp == "A"]$sd, 3); T[[paste0(p, "_sdLA")]] <- num(x[size == "l" & sp == "A"]$sd, 3)
  T[[paste0(p, "_sdSC")]] <- num(x[size == "s" & sp == "C"]$sd, 3); T[[paste0(p, "_sdLC")]] <- num(x[size == "l" & sp == "C"]$sd, 3)
  T[[paste0(p, "_fourA")]] <- num(x[size == "s" & sp == "A"]$sd / x[size == "l" & sp == "A"]$sd, 1)
  T[[paste0(p, "_fourB")]] <- num(x[size == "s" & sp == "B"]$sd / x[size == "l" & sp == "B"]$sd, 1)
  T[[paste0(p, "_fourC")]] <- num(x[size == "s" & sp == "C"]$sd / x[size == "l" & sp == "C"]$sd, 1)
}
# prior mean zero under selection
x <- g("select", "gsc0")
T$z_blo <- f(min(x$bias), 2); T$z_bhi <- f(max(x$bias), 2); T$z_covlo <- num(min(x$cov), 2); T$z_covhi <- num(max(x$cov), 2)
z <- function(sz, p) x[size == sz & sp == p]
T$z_n2se <- as.character(sum(abs(x$bias) > 2 * x$se_bias))
for (sz in c("s", "m", "l")) { y <- x[size == sz]; T[[paste0("z_mean_", sz)]] <- f(mean(y$bias), 2); T[[paste0("z_cov_", sz)]] <- num(mean(y$cov), 2) }
y <- g("none", "gsc0"); T$z0_maxb <- num(max(abs(y$bias)), 3); T$z0_covlo <- num(min(y$cov), 2); T$z0_covhi <- num(max(y$cov), 2)
# individual MSE
x <- g("select", "gsc1")
T$im_SA <- num(x[size == "s" & sp == "A"]$indiv, 2); T$im_LA <- num(x[size == "l" & sp == "A"]$indiv, 2); T$im_SC <- num(x[size == "s" & sp == "C"]$indiv, 2); T$im_LC <- num(x[size == "l" & sp == "C"]$indiv, 2)
T$s2_SA <- num(x[size == "s" & sp == "A"]$sigma2, 2); T$s2_MA <- num(x[size == "m" & sp == "A"]$sigma2, 2); T$s2_LA <- num(x[size == "l" & sp == "A"]$sigma2, 2)
y <- x[sp %in% c("B", "C")]; T$s2_lo <- num(min(y$sigma2), 2); T$s2_hi <- num(max(y$sigma2), 2)
# comparison at 40/100
cmp <- function(sc, m, p = "B", shp = NULL) { if (is.null(shp)) shp <- if (m %in% c("lme1", "lme2")) "linear" else "any"; s[scen == sc & size == "m" & sp == p & method == m & bin == "all" & shape == shp] }
sc_map <- c(none = "B", select = "B", growth = "B", outside = "B", wrongk = "B", visits = "B", dip = "B", dense = "C")
for (sc in names(sc_map)) for (m in c("gsc1", "gsc0", "did", "lme1", "lme2", "ife2", "ifecv", "mc")) {
  x <- cmp(sc, m, sc_map[[sc]])
  T[[paste0(m, "_", sc, "_b")]] <- f(x$bias, 2); T[[paste0(m, "_", sc, "_r")]] <- num(x$rmse, 2)
  T[[paste0(m, "_", sc, "_i")]] <- num(x$indiv, 2)
  if (x$ncov > 0) { T[[paste0(m, "_", sc, "_c")]] <- num(x$cov, 2); T[[paste0(m, "_", sc, "_w")]] <- num(x$width, 2) }
}
for (sc in names(sc_map)) { x <- cmp(sc, "ifecv", sc_map[[sc]]); T[[paste0("r0_", sc)]] <- as.character(round(100 * x$r0)); x <- cmp(sc, "ife2", sc_map[[sc]]); T[[paste0("kept_", sc)]] <- as.character(round(100 * x$kept)) }
for (p in c("A", "B", "C")) { x <- s[scen == "none" & size == "m" & sp == p & method == "ife2" & bin == "all"]; T[[paste0("kept", p)]] <- as.character(round(100 * x$kept)); T[[paste0("ifereps", p)]] <- as.character(x$reps); x <- s[scen == "none" & size == "m" & sp == p & method == "ifecv" & bin == "all"]; T[[paste0("cvreps", p)]] <- as.character(x$reps) }
x <- s[size == "m" & method %in% c("ife2", "ifecv") & sp == "A" & bin == "all"]; T$ifeA_lo <- as.character(min(x$reps)); T$ifeA_hi <- as.character(max(x$reps)); T$ifeA_R <- T$R_ife
x <- cmp("dense", "ifecv", "C"); T$cv_r_dense <- num(x$r_used, 1)
# outside the basis: width against rmse
x <- cmp("outside", "gsc1"); T$out_wmed <- num(x$width_med, 2)
for (p in c("A", "C")) for (sc in c("select", "growth", "outside")) for (m in c("gsc1", "gsc0", "did", "mc", "lme2")) { x <- cmp(sc, m, p); T[[paste0(m, "_", sc, "_", p, "_b")]] <- f(x$bias, 2); T[[paste0(m, "_", sc, "_", p, "_r")]] <- num(x$rmse, 2); if (x$ncov > 0) T[[paste0(m, "_", sc, "_", p, "_c")]] <- num(x$cov, 2) }
# LME effect shapes
for (sh in c("linear", "constant", "saturating", "transient")) for (m in c("lme1", "lme2")) { x <- cmp("none", m, "B", sh); T[[paste0(m, "_", sh, "_nb")]] <- f(x$bias, 2); x <- cmp("select", m, "B", sh); T[[paste0(m, "_", sh, "_sb")]] <- f(x$bias, 2); T[[paste0(m, "_", sh, "_sc")]] <- num(x$cov, 2) }
T$rshape <- as.character(cmp("none", "lme1", "B", "constant")$reps)
# LOOIC
for (sc in c("looic", "looic3")) {
  cc <- ch[scen == sc]; tot <- sum(cc$reps); p <- sc
  for (k in 1:3) { z <- cc[chosen_k == k]$reps; T[[paste0(p, "_ch", k)]] <- as.character(round(100 * (if (length(z)) z else 0) / tot)) }
  for (m in c("gsc1_k1", "gsc1_k2", "gsc1_k3", "gsc1_looic")) { x <- s[scen == sc & method == m & bin == "all"]; T[[paste0(p, "_", m, "_b")]] <- f(x$bias, 2); T[[paste0(p, "_", m, "_r")]] <- num(x$rmse, 2); T[[paste0(p, "_", m, "_c")]] <- num(x$cov, 2); T[[paste0(p, "_", m, "_w")]] <- num(x$width, 2) }
}
words <- c("none", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine")
T$n_2se_text <- { k <- as.integer(T$n_n2se); sprintf("%s of the nine estimates exceeded two Monte Carlo SEs", words[k + 1]) }
T$s_2se_text <- { k <- as.integer(T$s_n2se); sprintf("%s of the nine estimates exceeded two Monte Carlo SEs", words[k + 1]) }
T$ifeA_miss <- as.character(as.numeric(T$R_ife) - s[scen == "select" & size == "m" & method == "ife2" & sp == "A" & bin == "all"]$reps)
x01 <- s[scen %in% c("none", "select") & method == "gsc1" & bin == "0.1"]; T$t01_miss <- as.character(as.numeric(T$R) - min(x01$reps))
# the simpler-model rule of Section 2.3, on the first 100 datasets of each k-choice scenario
rc <- rbindlist(lapply(Sys.glob("rule_check_w*.csv"), fread))
pc <- function(cl, col, k) as.character(round(100 * mean(rc[cell == cl][[col]] == k)))
T$rule_n <- as.character(min(rc[, .N, by = cell]$N))
T$rule_21 <- pc(29, "z1", 2); T$rule_22 <- pc(29, "z2", 2); T$rule_11 <- pc(29, "z1", 1); T$rule_12 <- pc(29, "z2", 1)
T$rule_31 <- pc(30, "z1", 3); T$rule_32 <- pc(30, "z2", 3); T$rule_3k2_1 <- pc(30, "z1", 2); T$rule_3k2_2 <- pc(30, "z2", 2); T$rule_3k1_2 <- pc(30, "z2", 1)
print(rc[, .(n = .N, best2 = mean(best == 2), z1_1 = mean(z1 == 1), z1_2 = mean(z1 == 2), z1_3 = mean(z1 == 3), z2_1 = mean(z2 == 1), z2_2 = mean(z2 == 2), z2_3 = mean(z2 == 3)), by = cell])
T$rule_text <- if (min(rc[, .N, by = cell]$N) >= 100) sprintf("Under that rule, in the first %s datasets of each scenario, with true $k=2$ it chose $k=2$ in %s\\%% of datasets at one standard error and %s\\%% at two; with true $k=3$ it chose $k=3$ in %s\\%% and %s\\%%, and otherwise $k=2$, which has the bias shown in Table~\\ref{tab:looic}, so at two standard errors the rule can miss a component.", T$rule_n, T$rule_21, T$rule_22, T$rule_31, T$rule_32) else ""
write_json(T, "tokens_sim.json", auto_unbox = TRUE, pretty = TRUE)
cat(length(T), "tokens\n")
