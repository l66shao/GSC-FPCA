# Summarise the simulation. All cells are cut to the same number of replicates, the largest R that every cell reached
# (comparators keep their own caps), so the tables rest on one stated count. SIM_R overrides R.
library(data.table)
fs <- list.files("results2", pattern = "^w.*csv$", full.names = TRUE)
r <- rbindlist(lapply(fs, fread), fill = TRUE)
r <- unique(r, by = c("cell", "iter", "method", "shape", "bin"), fromLast = TRUE)
cells <- fread("cells2.csv")
# Every cell uses its first finished replicates, up to its target. A worker stopped mid-task can leave one replicate
# unfinished in a cell; that cell then uses the next replicate instead (stopping is unrelated to the outcome). Cells with
# the same target are cut to the same count, the largest that all of them reached, so each group of settings rests on
# one stated number. Comparators keep their own caps, which are by replicate number.
done <- r[, .(iter = sort(unique(iter))), by = cell]
done[, pos := seq_len(.N), by = cell]
cnt <- done[, .(k = .N), by = cell]
cc <- merge(cells[, .(cell, reps)], cnt, by = "cell", all.x = TRUE); cc[is.na(k), k := 0]
cc[, use := min(k, reps), by = cell]; cc[, use := min(use), by = reps]
r_all <- r
r <- merge(r, merge(done, cc[, .(cell, use)], by = "cell")[pos <= use, .(cell, iter)], by = c("cell", "iter"))
cat("Replicates used by target:\n"); print(unique(cc[, .(target = reps, used = use)]))
cat("Cells with a skipped replicate:", paste(unique(done[pos != iter]$cell), collapse = ","), "\n")
s <- r[, .(reps = sum(!is.na(err)), bias = mean(err, na.rm = TRUE), sd = sd(err, na.rm = TRUE), rmse = sqrt(mean(err^2, na.rm = TRUE)), cov = mean(lo <= 0 & hi >= 0, na.rm = TRUE), ncov = sum(!is.na(lo)), width = mean(hi - lo, na.rm = TRUE), width_med = median(hi - lo, na.rm = TRUE),
           indiv = mean(indiv_mse, na.rm = TRUE), indiv_sd = sd(indiv_mse, na.rm = TRUE), kept = mean(kept, na.rm = TRUE), n = mean(n, na.rm = TRUE), r_used = mean(r_used, na.rm = TRUE), r0 = mean(r_used == 0, na.rm = TRUE), sigma2 = mean(sigma2, na.rm = TRUE), nu1 = mean(nu1), nu2 = mean(nu2), secs = mean(secs, na.rm = TRUE), truth = mean(truth_linear, na.rm = TRUE)),
         by = .(cell, scen, size, sp, method, shape, bin)]
s[, se_bias := sd / sqrt(reps)]
s <- merge(s, cells[, .(cell, Ntr, Nco, lambda, n_min)], by = "cell")
fwrite(s, "summary2.csv")
# mixed models by effect shape, over the replicates that ran all four shapes
sh <- r[method %in% c("lme1", "lme2") & size == "m" & sp == "B" & scen %in% c("none", "select") & bin == "all" & iter <= 100, .(reps = sum(!is.na(err)), bias = mean(err, na.rm = TRUE), rmse = sqrt(mean(err^2, na.rm = TRUE)), cov = mean(lo <= 0 & hi >= 0, na.rm = TRUE)), by = .(scen, method, shape)]
fwrite(sh, "shape100.csv")
ch <- r[method == "gsc1_looic" & bin == "all", .(reps = .N), by = .(cell, scen, chosen_k)]
fwrite(ch, "looic_choice.csv")
cat("tasks done:", nrow(unique(r_all[, .(cell, iter)])), " rows", nrow(r_all), "\n")
print(merge(cells[, .(cell, scen, size, sp, reps)], done[, .(finished = .N), by = cell], by = "cell", all.x = TRUE)[, .(cell, scen, size, sp, target = reps, finished)], nrows = 40)
