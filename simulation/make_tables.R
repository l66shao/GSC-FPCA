# Writes LaTeX tables for Section 3 and Appendix D from summary2.csv, plus a CSV of numbers used in the prose.
library(data.table)
s <- fread("summary2.csv")
ch <- fread("looic_choice.csv")
neg <- function(x, d = 3) { v <- sprintf(paste0("%.", d, "f"), abs(x)); ifelse(is.na(x), "--", ifelse(x < 0 & as.numeric(v) != 0, paste0("$-$", v), v)) }
cov3 <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
num <- function(x, d = 2) ifelse(is.na(x), "--", sprintf(paste0("%.", d, "f"), x))
szlab <- c(s = "20/50", m = "40/100", l = "80/200"); splab <- c(A = "2 / 2.5", B = "3 / 4.3", C = "4 / 6.9")
get <- function(sc, sz, p, m, b = "all", shp = "any") { x <- s[scen == sc & size == sz & sp == p & method == m & bin == b & shape == shp]; if (nrow(x)) x else data.table(bias = NA_real_, sd = NA_real_, cov = NA_real_, indiv = NA_real_, width = NA_real_, rmse = NA_real_, kept = NA_real_, reps = 0, se_bias = NA_real_, indiv_sd = NA_real_, n = NA_real_, r_used = NA_real_) }
# ---- Table 1 ----
rows <- character()
for (sz in c("s", "m", "l")) {
  for (p in c("A", "B", "C")) {
    a <- get("none", sz, p, "gsc1"); b <- get("select", sz, p, "gsc1"); z <- get("select", sz, p, "gsc0")
    rows <- c(rows, sprintf("%s & %s & %s\\,(%s) & %s & %s\\,(%s) & %s & %s & %s & %s\\,(%s) \\\\", szlab[sz], splab[p], neg(a$bias), num(a$sd, 3), cov3(a$cov), neg(b$bias), num(b$sd, 3), cov3(b$cov), neg(z$bias), cov3(z$cov), num(b$indiv, 3), num(b$indiv_sd, 3)))
  }
  if (sz != "l") rows <- c(rows, "\\midrule")
}
writeLines(rows, "tab_main_rows.tex")
# ---- Table at scaled time 0.1 (Appendix D) ----
rows <- character()
for (sz in c("s", "m", "l")) {
  for (p in c("A", "B", "C")) {
    a <- get("none", sz, p, "gsc1", "0.1"); b <- get("select", sz, p, "gsc1", "0.1"); z <- get("select", sz, p, "gsc0", "0.1")
    rows <- c(rows, sprintf("%s & %s & %s\\,(%s) & %s & %s\\,(%s) & %s & %s & %s \\\\", szlab[sz], splab[p], neg(a$bias), num(a$sd, 3), cov3(a$cov), neg(b$bias), num(b$sd, 3), cov3(b$cov), neg(z$bias), cov3(z$cov)))
  }
  if (sz != "l") rows <- c(rows, "\\midrule")
}
writeLines(rows, "tab_t01_rows.tex")
# ---- Table 2: comparison at 40/100, sparsity B, in two panels so that it fits the text width ----
scens <- list(c("none", "B"), c("select", "B"), c("growth", "B"), c("outside", "B"), c("wrongk", "B"), c("visits", "B"), c("dip", "B"), c("dense", "C"))
meths <- list(gsc1 = "GSC-FPCA", gsc0 = "GSC-FPCA, mean 0", did = "DiD", lme1 = "LME-1", lme2 = "LME-2", ife2 = "IFE, $r=2$", ifecv = "IFE, CV", mc = "MC")
row_for <- function(m, f, sc = scens) { shp <- if (m %in% c("lme1", "lme2")) "linear" else "any"; paste(sapply(sc, function(a) { x <- get(a[1], "m", a[2], m, "all", shp); f(x) }), collapse = " & ") }
panel <- function(sc, cov_m) {
  A <- sapply(names(meths), function(m) sprintf("%s & %s \\\\", meths[[m]], row_for(m, function(x) sprintf("%s\\,(%s)", neg(x$bias, 2), num(x$rmse, 2)), sc)))
  B <- sapply(cov_m, function(m) sprintf("%s & %s \\\\", meths[[m]], row_for(m, function(x) if (is.na(x$cov) || x$ncov == 0) "--" else sprintf("%s\\,(%s)", num(x$cov, 2), num(x$width, 2)), sc)))
  c("\\multicolumn{5}{l}{\\emph{Bias (RMSE)}} \\\\", A, "\\addlinespace", "\\multicolumn{5}{l}{\\emph{Coverage of nominal 95\\% interval (average width)}} \\\\", B)
}
writeLines(panel(scens[1:4], c("gsc1", "gsc0", "did", "lme1", "lme2", "mc")), "tab_compare_a_rows.tex")
writeLines(panel(scens[5:8], c("gsc1", "gsc0", "did", "lme1", "lme2")), "tab_compare_b_rows.tex")
# individual-level MSE and units kept by IFE (Appendix D)
C <- sapply(c("gsc1", "gsc0", "did", "ife2", "ifecv", "mc"), function(m) sprintf("%s & %s \\\\", meths[[m]], row_for(m, function(x) num(x$indiv, 2))))
K <- sprintf("Share kept by IFE & %s \\\\", row_for("ife2", function(x) num(x$kept, 2)))
writeLines(c(C, "\\midrule", K), "tab_compare_indiv_rows.tex")
# ---- comparison at sparsity A and C (Appendix D) ----
rows <- character()
scens2 <- c("none", "select", "growth", "outside")
for (p in c("A", "C")) {
  rows <- c(rows, sprintf("\\multicolumn{5}{l}{\\emph{Sparsity %s (minimum/average $n_i$: %s)}} \\\\", p, splab[p]))
  for (m in names(meths)) { shp <- if (m %in% c("lme1", "lme2")) "linear" else "any"
    rows <- c(rows, sprintf("%s & %s \\\\", meths[[m]], paste(sapply(scens2, function(sc) { x <- get(sc, "m", p, m, "all", shp); sprintf("%s\\,(%s)", neg(x$bias, 2), num(x$rmse, 2)) }), collapse = " & ")))
  }
  if (p == "A") rows <- c(rows, "\\midrule")
}
writeLines(rows, "tab_compare_AC_rows.tex")
# ---- LME by effect shape, over the same datasets for all four shapes ----
s100 <- fread("shape100.csv")
rows <- character()
for (sh in c("linear", "constant", "saturating", "transient")) {
  v <- sapply(c("none", "select"), function(sc) sapply(c("lme1", "lme2"), function(m) { x <- s100[scen == sc & method == m & shape == sh]; sprintf("%s\\,(%s) & %s", neg(x$bias, 2), num(x$rmse, 2), num(x$cov, 2)) }))
  rows <- c(rows, sprintf("%s & %s \\\\", tools::toTitleCase(sh), paste(as.vector(v), collapse = " & ")))
}
writeLines(rows, "tab_lme_rows.tex")
# ---- k by LOOIC ----
rows <- character()
for (sc in c("looic", "looic3")) {
  tk <- if (sc == "looic") 2 else 3
  cc <- ch[scen == sc]; tot <- sum(cc$reps); pr <- sapply(1:3, function(k) { z <- cc[chosen_k == k]$reps; if (length(z)) z / tot else 0 })
  for (m in c("gsc1_k1", "gsc1_k2", "gsc1_k3", "gsc1_looic")) {
    x <- s[scen == sc & method == m & bin == "all"]
    lab <- if (m == "gsc1_looic") "$k$ by LOOIC" else paste0("$k=", substr(m, 7, 7), "$")
    rows <- c(rows, sprintf("%s & %s & %s & %s\\,(%s) & %s & %s & %s \\\\", tk, lab, if (m == "gsc1_looic") "--" else num(pr[as.integer(substr(m, 7, 7))], 2), neg(x$bias, 3), num(x$rmse, 3), num(x$cov, 3), num(x$width, 2), num(x$indiv, 2)))
  }
  if (sc == "looic") rows <- c(rows, "\\midrule")
}
writeLines(rows, "tab_looic_rows.tex")
cat("reps by scenario:\n"); print(s[bin == "all" & method == "gsc1", .(reps = max(reps)), by = scen]); print(s[bin == "all" & method %in% c("lme1", "ife2", "ifecv", "mc") & size == "m" & sp == "B" & scen %in% c("none", "select") & shape %in% c("any", "linear"), .(scen, method, reps, ncov)])
