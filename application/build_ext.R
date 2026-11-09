# Rebuild the analysis rows from FullData_Y9 with the steps of extract_ncanda_data_25-10-17.R, confirm they match
# ncanda_dat.csv, and add ROI outcomes and richer covariates. Prints aggregate checks only.
library(data.table)
F <- Sys.getenv("NCANDA_FULLDATA", "data/FullData_Y9.csv")
hd <- names(fread(F, nrows = 0))
gv <- grep("_grayvol\\.[xy]$", hd, value = TRUE)
keep <- c("subject", "arm", "visit", "visit_age", "site", "sex", "race", "hispanic", "exceeds_bl_drinking", "sri24_icv_brain_volume", "sri24_icv_brain_gm", "sri24_icv_brain_wm",
  "youthreport1_cddr27", "youthreport1_cddr30", "youthreport1_cddr4a", "youthreport1_cddr4b", "youthreport1_cddr55a", "youthreport1_cddr55b", "youthreport1_cddr19", "youthreport1_cddr20",
  "ses_parent_yoe", "ses_parent_faminc", "ses_parent_score", "parent_edu", "fh_alc_density", "fh_drug_density", "pds_score", "pds_pubcat", "bmi_zscore", "family_id", "scanner", "scanner_model", "mri_structural_age",
  "ssaga_youth_externalizing", "highrisk_ssaga_extern_tot", "cesd_score", gv)
keep <- intersect(keep, hd)
d <- fread(F, select = keep)
cat("rows", nrow(d), "subjects", uniqueN(d$subject), "grayvol cols", length(gv), "\n")
d <- as.data.frame(d)
d$binge_YN <- d$youthreport1_cddr27; d$binge <- d$youthreport1_cddr30
d$binge[d$binge_YN == 0 & !is.na(d$binge_YN)] <- 0
d <- d[d$exceeds_bl_drinking == "N", ]
cat("after exceeds filter: subjects", length(unique(d$subject)), "\n")
d$binge[is.na(d$binge & d$visit == "baseline")] <- 0
d <- d[!is.na(d$binge), ]
cat("after NA binge filter: rows", nrow(d), "subjects", length(unique(d$subject)), "\n")
vm <- c(baseline = 1, setNames(2:10, paste0("followup_", 1:9, "y")))
d$visit <- vm[d$visit]
d$id <- match(d$subject, unique(d$subject))
# gray matter: left/right average, as in the extraction script
rois <- sub("\\.x$", "", grep("_grayvol\\.x$", gv, value = TRUE))
cat("n ROI columns x:", length(rois), "\n")
# the extraction script used grayvol columns 3:36 paired with +34; check pairing by name
pair_ok <- all(rois %in% sub("\\.y$", "", grep("_grayvol\\.y$", gv, value = TRUE)))
cat("pairs ok", pair_ok, "\n"); print(rois)
for (r in rois) d[[r]] <- (d[[paste0(r, ".x")]] + d[[paste0(r, ".y")]]) / 2
d$sf_left <- d[["superiorfrontal_grayvol.x"]]; d$sf_right <- d[["superiorfrontal_grayvol.y"]]
d[paste0(rois, ".x")] <- NULL; d[paste0(rois, ".y")] <- NULL
o <- read.csv(Sys.getenv("NCANDA_DAT", "data/ncanda_dat.csv"))
cat("ncanda_dat rows", nrow(o), "ids", length(unique(o$id)), "\n")
m <- merge(o, d, by = c("id", "visit"), suffixes = c("", ".new"), all.x = TRUE)
cat("rows matched", sum(!is.na(m$subject)), "of", nrow(o), "\n")
chk <- function(a, b) c(n = sum(!is.na(a) & !is.na(b)), maxdiff = max(abs(a - b), na.rm = TRUE), onlyA = sum(!is.na(a) & is.na(b)), onlyB = sum(is.na(a) & !is.na(b)))
print(rbind(age = chk(m$visit_age, m$visit_age.new), icv = chk(m$sri24_icv_brain_volume, m$sri24_icv_brain_volume.new), sf = chk(m$superiorfrontal_grayvol, m$superiorfrontal_grayvol.new), binge = chk(m$binge, m$binge.new)))
print(colSums(is.na(m[, intersect(c("ses_parent_yoe", "ses_parent_faminc", "ses_parent_score", "parent_edu", "fh_alc_density", "fh_drug_density", "pds_score", "bmi_zscore", "youthreport1_cddr4a", "youthreport1_cddr4b", "youthreport1_cddr55a", "youthreport1_cddr55b", "scanner", "family_id"), names(m))])))
saveRDS(m, "ext_all.rds"); saveRDS(rois, "rois.rds")
