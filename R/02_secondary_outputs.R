# Secondary analyses and supplementary tables (S6 / S7b / S8 / S8b / S10 / S12 / S13)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(pROC)
})

# Auto-detect PROJECT_ROOT relative to this script's location
find_project_root <- function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  file_arg <- cargs[grep("^--file=", cargs)]
  if (length(file_arg) == 1L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
    script_dir  <- dirname(script_path)
  } else {
    script_dir <- getwd()
  }
  root <- script_dir
  for (i in 1:6) {
    if (file.exists(file.path(root, "R", "00_functions.R"))) {
      return(normalizePath(root, winslash = "/"))
    }
    parent <- dirname(root)
    if (parent == root) break   # reached filesystem root
    root <- parent
  }
  stop("Cannot locate project root (R/00_functions.R not found).")
}
PROJECT_ROOT <- find_project_root()
source(file.path(PROJECT_ROOT, "R", "00_functions.R"))

TABLES <- file.path(PROJECT_ROOT, "output", "tables")
dir.create(TABLES, showWarnings = FALSE, recursive = TRUE)

dat <- fread(file.path(PROJECT_ROOT, "data", "feature_matrix_with_clusters.csv"))

# Time-to-event: 28-day death, hospital LOS censored at 28d
parse_los <- function(x) {
  parts <- strsplit(x, " ")[[1]]
  days <- as.numeric(parts[1])
  hms  <- as.numeric(strsplit(parts[3], ":")[[1]])
  days + hms[1]/24 + hms[2]/1440 + hms[3]/86400
}
los_days <- sapply(dat$hosp_los_interval, parse_los)
dat[, time := pmin(los_days, 28)]
dat[, event := death_28day]
dat[, wbc_high := as.integer(wbc_baseline > 15.3)]

vars <- c("admission_age", "sex_f", "sofa_sepsis3", "charlson",
          "has_vent", "has_vaso", "has_crrt",
          "immunocompromised",
          "has_glycopeptide", "has_carbapenem", "has_cephalosporin")
cols <- c("time", "event", "km_class", "wbc_baseline", "wbc_high", vars)
dd <- na.omit(dat[, ..cols])
cat(sprintf("Complete cases: %d (%.1f%% of %d)\n", nrow(dd), nrow(dd)/nrow(dat)*100, nrow(dat)))

f_km   <- as.formula(paste("Surv(time, event) ~ km_class +",     paste(vars, collapse=" + ")))
f_wbc  <- as.formula(paste("Surv(time, event) ~ wbc_baseline +", paste(vars, collapse=" + ")))
f_hi   <- as.formula(paste("Surv(time, event) ~ wbc_high +",     paste(vars, collapse=" + ")))
f_both <- as.formula(paste("Surv(time, event) ~ wbc_baseline + km_class +", paste(vars, collapse=" + ")))
f_hi2  <- as.formula(paste("Surv(time, event) ~ wbc_high + km_class +",     paste(vars, collapse=" + ")))

m_km  <- coxph(f_km,  data = dd)
m_wbc <- coxph(f_wbc, data = dd)
m_hi  <- coxph(f_hi,  data = dd)
m_both <- coxph(f_both, data = dd)
m_hi2  <- coxph(f_hi2, data = dd)

c_km  <- summary(m_km)$concordance[1]
c_wbc <- summary(m_wbc)$concordance[1]
c_hi  <- summary(m_hi)$concordance[1]

cat(sprintf("\n── C-index ──────────────────────────────\n"))
cat(sprintf("Model A (km_class):      %.4f\n", c_km))
cat(sprintf("Model B (wbc_baseline):   %.4f  (Δ=%+.4f vs A)\n", c_wbc, c_wbc-c_km))
cat(sprintf("Model C (wbc>15.3):       %.4f  (Δ=%+.4f vs A)\n", c_hi,  c_hi-c_km))

# LRT: incremental value of km_class
lrt_wbc <- anova(m_wbc, m_both, test = "LRT")
lrt_hi  <- anova(m_hi,  m_hi2,  test = "LRT")
p_wbc <- lrt_wbc[["Pr(>|Chi|)"]][2]
p_hi  <- lrt_hi[["Pr(>|Chi|)"]][2]

cat(sprintf("\n── LRT (km_class adds to static WBC) ────\n"))
cat(sprintf("+ km_class over wbc_baseline: χ²=%.2f, df=%d, p=%.2e\n",
            lrt_wbc$Chisq[2], lrt_wbc$Df[2], p_wbc))
cat(sprintf("+ km_class over wbc>15.3:     χ²=%.2f, df=%d, p=%.2e\n",
            lrt_hi$Chisq[2],  lrt_hi$Df[2],  p_hi))

res <- data.table(
  Model            = c("km_class","wbc_baseline","wbc>15.3"),
  C_index          = c(c_km, c_wbc, c_hi),
  Delta_vs_km      = c(0, c_wbc-c_km, c_hi-c_km),
  LRT_add_km_chisq = c(NA, lrt_wbc$Chisq[2], lrt_hi$Chisq[2]),
  LRT_add_km_df    = c(NA, lrt_wbc$Df[2],   lrt_hi$Df[2]),
  LRT_add_km_p     = c(NA, p_wbc, p_hi)
)
# LRT_add_km_* lives on the "wbc_baseline" row (the "wbc>15.3" row is the same test)
fwrite(res, file.path(TABLES, "tableS8_trajectory_vs_static_wbc.csv"))
cat("\nSaved.\n")


if (!file.exists(file.path(PROJECT_ROOT, "data", "feature_matrix.csv")))
  stop("data/feature_matrix.csv not found - run the SQL extraction first")
cat("PROJECT_ROOT:", PROJECT_ROOT, "\n")
cat("R:", R.version.string, "\n")

operating_stats <- function(y1, wbc, cutoff = 15.3) {
  # y1: 1 = C1 (positive); wbc: numeric baseline WBC; fixed rule: wbc > cutoff
  ok <- !is.na(y1) & !is.na(wbc)
  y1 <- as.integer(y1[ok]); wbc <- wbc[ok]
  pred <- as.integer(wbc > cutoff)
  tp <- sum(pred == 1 & y1 == 1); fp <- sum(pred == 1 & y1 == 0)
  tn <- sum(pred == 0 & y1 == 0); fn <- sum(pred == 0 & y1 == 1)
  sens <- tp / max(1, tp + fn); spec <- tn / max(1, tn + fp)
  ppv  <- tp / max(1, tp + fp); npv  <- tn / max(1, tn + fn)
  data.table(n = length(y1), prev_pct = mean(y1) * 100,
             sens_pct = sens * 100, spec_pct = spec * 100,
             ppv_pct = ppv * 100, npv_pct = npv * 100,
             wbc_gt_n = sum(pred), wbc_le_n = sum(pred == 0))
}

auc_youden <- function(y1, wbc) {
  ok <- !is.na(y1) & !is.na(wbc)
  r <- roc(as.integer(y1[ok]), wbc[ok], quiet = TRUE)
  auc_v <- as.numeric(auc(r))
  youden <- r$sensitivities + r$specificities - 1
  oi <- which.max(youden)
  data.table(auc = auc_v, youden_cutoff = r$thresholds[oi],
             youden_sens = r$sensitivities[oi] * 100,
             youden_spec = r$specificities[oi] * 100)
}

hr_direct <- function(fit, term) {
  b <- coef(fit)[term]; ci <- confint(fit)[term, ]
  c(hr = unname(exp(b)), lo = unname(exp(ci[1])), hi = unname(exp(ci[2])))
}

# --- B1: MIMIC reference ---
cat("\n========== B1: MIMIC ==========\n")
f <- fread(file.path(PROJECT_ROOT, "data", "feature_matrix_with_clusters.csv"))
cat("MIMIC n =", nrow(f), "\n")
st_mimic <- operating_stats(f$km_class == 1, f$wbc_baseline, 15.3)
au_mimic <- auc_youden(f$km_class == 1, f$wbc_baseline)
st_mimic[, cutoff := 15.3][, cohort := "MIMIC"]
cat("AUC:", au_mimic$auc, " youden cutoff:", au_mimic$youden_cutoff, "\n")
print(st_mimic)
# mortality by threshold (rule-in / rule-out) - 28-day
f[, wbc_hi := as.integer(wbc_baseline > 15.3)]
mort_by <- f[!is.na(wbc_hi), .(n = .N, mort_28d_pct = mean(death_28day) * 100,
  c1_pct = mean(km_class == 1) * 100), by = wbc_hi][order(-wbc_hi)]
print(mort_by)

# --- B1: eICU ---
cat("\n========== B1: eICU ==========\n")
eicu_l <- fread(file.path(PROJECT_ROOT, "data", "eicu_blood_long.csv"))
eicu_w <- build_wbc_windows(eicu_l, "patientunitstayid")
eicu_a <- fread(file.path(PROJECT_ROOT, "data", "eicu_cluster_assignment.csv"))
eicu <- merge(eicu_a, eicu_w, by = "patientunitstayid")
cat("eICU with assignment: n =", nrow(eicu),
    " baseline WBC missing:", sum(is.na(eicu$wbc_baseline)), "\n")
# de novo labels (Method A) primary
st_eicu_a <- operating_stats(eicu$km_class_recluster == 1, eicu$wbc_baseline, 15.3)
au_eicu_a <- auc_youden(eicu$km_class_recluster == 1, eicu$wbc_baseline)
st_eicu_a[, cutoff := 15.3][, cohort := "eICU (de novo labels)"]
cat("AUC(de novo):", au_eicu_a$auc, "\n"); print(st_eicu_a)
# Method B (MIMIC centroid) labels
st_eicu_b <- operating_stats(eicu$km_class_mimic == 1, eicu$wbc_baseline, 15.3)
au_eicu_b <- auc_youden(eicu$km_class_mimic == 1, eicu$wbc_baseline)
st_eicu_b[, cutoff := 15.3][, cohort := "eICU (MIMIC-centroid labels)"]
cat("AUC(mimic labels):", au_eicu_b$auc, "\n"); print(st_eicu_b)
# mortality by threshold (hospital death)
eicu[, wbc_hi := as.integer(wbc_baseline > 15.3)]
mort_eicu <- eicu[!is.na(wbc_hi) & !is.na(hosp_death),
  .(n = .N, hosp_mort_pct = mean(hosp_death) * 100,
    c1_pct = mean(km_class_recluster == 1) * 100), by = wbc_hi][order(-wbc_hi)]
print(mort_eicu)

# --- B1: NWICU (de novo clusters, same rules as R/01) ---
cat("\n========== B1: NWICU ==========\n")
nw_l  <- fread(file.path(PROJECT_ROOT, "data", "nwicu_blood_long.csv"))
nw_c  <- fread(file.path(PROJECT_ROOT, "data", "nwicu_cohort.csv"))
nw_l[, `:=`(stay_id = as.integer(stay_id), hrs = as.numeric(hrs), wbc = as.numeric(wbc))]
nw_f  <- build_wbc_windows(nw_l, "stay_id")
X <- as.matrix(nw_f[, ..CLUSTER_FEATURES]); mode(X) <- "numeric"
for (j in seq_len(ncol(X))) {
  na_idx <- is.na(X[, j]); if (any(na_idx)) X[na_idx, j] <- median(X[, j], na.rm = TRUE)
}
X_scaled <- scale(X); X_scaled[!is.finite(X_scaled)] <- 0
set.seed(2026)
kf <- fit_kmeans(X_scaled, 2)
al <- align_labels_by_wbc(kf$km$cluster, X[, "wbc_baseline"])
nw_f[, km_class_recluster := al$labels]
model_path <- file.path(PROJECT_ROOT, "output", "model", "mimic_cluster_model.rds")
if (!file.exists(model_path))
  stop("MIMIC cluster model not found: ", model_path,
       "\nThis .rds is produced by R/01_cohort_analyses.R; run that script first.")
mm <- readRDS(model_path)
nw_f[, km_class_mimic := predict_cluster(nw_f, mm)]
nw <- merge(nw_c, nw_f[, .(stay_id, km_class_recluster, km_class_mimic)], by = "stay_id")
cat("NWICU n =", nrow(nw), "\n")
st_nw_a <- operating_stats(nw$km_class_recluster == 1, nw$wbc_baseline, 15.3)
au_nw_a <- auc_youden(nw$km_class_recluster == 1, nw$wbc_baseline)
st_nw_a[, cutoff := 15.3][, cohort := "NWICU (de novo labels)"]
cat("AUC(de novo):", au_nw_a$auc, "\n"); print(st_nw_a)
st_nw_b <- operating_stats(nw$km_class_mimic == 1, nw$wbc_baseline, 15.3)
au_nw_b <- auc_youden(nw$km_class_mimic == 1, nw$wbc_baseline)
st_nw_b[, cutoff := 15.3][, cohort := "NWICU (MIMIC-centroid labels)"]
cat("AUC(mimic labels):", au_nw_b$auc, "\n"); print(st_nw_b)
# mortality by threshold (28-day)
nw[, wbc_hi := as.integer(wbc_baseline > 15.3)]
mort_nw <- nw[!is.na(wbc_hi) & !is.na(death_28day),
  .(n = .N, mort_28d_pct = mean(death_28day) * 100,
    c1_pct = mean(km_class_recluster == 1) * 100), by = wbc_hi][order(-wbc_hi)]
print(mort_nw)

b1_out <- rbindlist(list(st_mimic, st_eicu_a, st_eicu_b, st_nw_a, st_nw_b),
                    fill = TRUE)
b1_out[, auc := c(au_mimic$auc, au_eicu_a$auc, au_eicu_b$auc,
                  au_nw_a$auc, au_nw_b$auc)]
fwrite(b1_out, file.path(TABLES, "_b1_threshold_validation.csv"))
mort_comb <- rbindlist(list(
  mort_by[, .(cohort = "MIMIC", wbc_hi, n, mortality_pct = mort_28d_pct, c1_pct)],
  mort_eicu[, .(cohort = "eICU", wbc_hi, n, mortality_pct = hosp_mort_pct, c1_pct)],
  mort_nw[, .(cohort = "NWICU", wbc_hi, n, mortality_pct = mort_28d_pct, c1_pct)]))
fwrite(mort_comb, file.path(TABLES, "_b1_threshold_mortality.csv"))

# --- B2: RERI in TRUE culture-negative (>=1 blood culture, no pathogen) ---
cat("\n========== B2: true culture-negative RERI ==========\n")
cat("cols present:", "n_blood_cultures" %in% names(f), " is_cn:", "is_cn" %in% names(f), "\n")
if ("n_blood_cultures" %in% names(f) && "is_cn" %in% names(f)) {
  tc <- f[is_cn == TRUE & n_blood_cultures >= 1 & !is.na(lactate_high)]
  cat("True culture-negative with lactate: n =", nrow(tc), "\n")
  tc2 <- tc[!is.na(death_28day)]
  cat("  C1 n =", sum(tc2$km_class == 1), " C2 n =", sum(tc2$km_class == 2), "\n")
  cat("  C1 mortality:", mean(tc2$death_28day[tc2$km_class == 1]) * 100,
      " C2 mortality:", mean(tc2$death_28day[tc2$km_class == 2]) * 100, "\n")
  tc2[, `:=`(A = as.integer(km_class == 1), B = as.integer(lactate_high))]
  r_crude <- compute_reri(tc2, "death_28day", "A", "B")
  r_adj   <- compute_reri(tc2, "death_28day", "A", "B",
             covars = c("admission_age", "sex_f", "sofa_sepsis3", "charlson"))
  r_out <- data.table(Measure = r_crude$Measure, Crude = r_crude$Value,
                      Adjusted = r_adj$Value)
  r_out[, n := nrow(tc2)]
  print(r_out)
  fwrite(r_out, file.path(TABLES, "_b2_true_cn_reri.csv"))
  cell <- tc2[, .(n = .N, mort_pct = mean(death_28day) * 100),
              by = .(C1 = km_class == 1, Lac4 = lactate_high)][order(-C1, -Lac4)]
  print(cell)
  fwrite(cell, file.path(TABLES, "_b2_true_cn_cells.csv"))
} else {
  cat("B2 skipped: required columns missing\n")
}

cat("\n========== ALL DONE ==========\n")

# ---------------- B3 ----------------
cat("\n========== B3: level vs shape (MIMIC Cox) ==========\n")
f <- fread(file.path(PROJECT_ROOT, "data", "feature_matrix_with_clusters.csv"))
parse_los <- function(x) {
  parts <- strsplit(x, " ")[[1]]
  days <- as.numeric(parts[1])
  hms  <- as.numeric(strsplit(parts[3], ":")[[1]])
  days + hms[1] / 24 + hms[2] / 1440 + hms[3] / 86400
}
f[, t_ev := pmin(sapply(hosp_los_interval, parse_los), 28)]
f[, ev := death_28day]
vars <- c("admission_age", "sex_f", "sofa_sepsis3", "charlson",
          "has_vent", "has_vaso", "has_crrt", "immunocompromised",
          "has_glycopeptide", "has_carbapenem", "has_cephalosporin")
cols <- c("t_ev", "ev", "km_class", "wbc_baseline", vars)
dd <- na.omit(f[, ..cols])
cat("Complete cases:", nrow(dd), "/", nrow(f), "\n")
dd[, c1 := as.integer(km_class == 1)]
cv <- paste(vars, collapse = " + ")
m1 <- coxph(as.formula(sprintf("Surv(t_ev, ev) ~ c1 + %s", cv)), data = dd)
m2 <- coxph(as.formula(sprintf("Surv(t_ev, ev) ~ c1 + wbc_baseline + %s", cv)), data = dd)
m3 <- coxph(as.formula(sprintf("Surv(t_ev, ev) ~ wbc_baseline + %s", cv)), data = dd)
hr_direct <- function(fit, term) {
  b <- coef(fit)[term]; ci <- confint(fit)[term, ]
  c(hr = unname(exp(b)), lo = unname(exp(ci[1])), hi = unname(exp(ci[2])))
}
r1 <- hr_direct(m1, "c1"); r2 <- hr_direct(m2, "c1")
lrt12 <- anova(m1, m2, test = "LRT")
b3 <- data.table(
  model = c("km_class (trajectory), no baseline WBC",
            "km_class + baseline WBC",
            "baseline WBC only (km_class omitted)"),
  hr_c1_vs_c2 = c(r1[["hr"]], r2[["hr"]], NA_real_),
  ci_lo = c(r1[["lo"]], r2[["lo"]], NA_real_),
  ci_hi = c(r1[["hi"]], r2[["hi"]], NA_real_),
  c_index = c(summary(m1)$concordance[1], summary(m2)$concordance[1],
              summary(m3)$concordance[1]),
  n = nrow(dd))
cat("HR km_class (no baseline):", round(r1[["hr"]], 3),
    "(", round(r1[["lo"]], 3), "-", round(r1[["hi"]], 3), ")\n")
cat("HR km_class (with baseline WBC):", round(r2[["hr"]], 3),
    "(", round(r2[["lo"]], 3), "-", round(r2[["hi"]], 3), ")\n")
cat("LRT adding baseline WBC to km_class model: chisq =",
    round(lrt12$Chisq[2], 2), " p =", format.pval(lrt12[["Pr(>|Chi|)"]][2]), "\n")
sl_cols <- c("km_class", "wbc_baseline", "wbc_t24h", "wbc_t48h", "wbc_t72h")
sl <- na.omit(f[, ..sl_cols])
sl[, slope := (wbc_t72h - wbc_baseline) / 3]
b3_shape <- sl[, .(n = .N, slope_mean = mean(slope), slope_med = median(slope),
                   base_mean = mean(wbc_baseline), base_med = median(wbc_baseline)),
               by = km_class][order(km_class)]
setnames(b3_shape, "km_class", "c1")
print(b3); print(b3_shape)
fwrite(b3, file.path(TABLES, "_b3_level_vs_shape.csv"))
fwrite(b3, file.path(TABLES, "tableS8b_level_vs_shape.csv"))
fwrite(b3_shape, file.path(TABLES, "_b3_slope_by_c1.csv"))

# ---------------- cluster stability ----------------
cat("\n========== ST: cluster stability refresh ==========\n")

# --- (a) Primary analysis matrix ---
mice_res <- impute_wbc_mice(f, features = CLUSTER_FEATURES,
                            aux = "n_tp_wbc", m = 5, seed = 2026)
Xmice <- mice_res$per_imp[[1]]
sc_mice <- scale(Xmice); Xmice_s <- sc_mice; Xmice_s[!is.finite(Xmice_s)] <- 0
ref <- f$km_class
kf <- fit_kmeans(Xmice_s, 2)                      # nstart = 50, seed = 2026
ari_repro <- adjustedRandIndex(kf$km$cluster, ref)
cat(sprintf("primary-matrix reproduction check: ARI vs stored km_class = %.6f\n",
            ari_repro))

# --- (b) Initialization / seed sensitivity (MICE primary base) ---
set.seed(2026)
init_aris <- vapply(1:25, function(i) {
  cl <- kmeans(Xmice_s, centers = 2, nstart = 1, iter.max = 100)$cluster
  adjustedRandIndex(cl, kf$km$cluster)
}, numeric(1))
cat(sprintf("initialization stability (nstart=1, 25 re-inits): mean %.3f range %.3f-%.3f\n",
            mean(init_aris), min(init_aris), max(init_aris)))

# --- (c) Bootstrap: MICE base (subsampling robustness) ---
set.seed(2026)
nboot <- 500
boot_mice <- numeric(nboot)
for (b in seq_len(nboot)) {
  idx <- sort(sample(nrow(Xmice_s), round(nrow(Xmice_s) * 0.7)))
  km <- kmeans(Xmice_s[idx, ], centers = 2, nstart = 25, iter.max = 100)
  boot_mice[b] <- adjustedRandIndex(km$cluster, kf$km$cluster[idx])
}
cat(sprintf("bootstrap (MICE base): mean %.3f  95%% %.3f-%.3f\n",
            mean(boot_mice), quantile(boot_mice, 0.025), quantile(boot_mice, 0.975)))

# --- (d) Sensitivity: median-imputation base (imputation + subsampling combined) ---
Xm <- as.matrix(f[, ..CLUSTER_FEATURES]); mode(Xm) <- "numeric"
for (j in seq_len(ncol(Xm))) {
  na_idx <- is.na(Xm[, j]); if (any(na_idx)) Xm[na_idx, j] <- median(Xm[, j], na.rm = TRUE)
}
Xm_s <- scale(Xm); Xm_s[!is.finite(Xm_s)] <- 0
set.seed(2026)
med_aris <- vapply(1:100, function(i) {
  km <- kmeans(Xm_s, centers = 2, nstart = 50, iter.max = 100)
  lab <- align_labels_by_wbc(km$cluster, Xm[, 1])$labels
  adjustedRandIndex(lab, ref)
}, numeric(1))
set.seed(2026)
boot_med <- numeric(nboot)
for (b in seq_len(nboot)) {
  idx <- sort(sample(nrow(Xm_s), round(nrow(Xm_s) * 0.7)))
  km <- kmeans(Xm_s[idx, ], centers = 2, nstart = 25, iter.max = 100)
  lab <- align_labels_by_wbc(km$cluster, Xm[idx, 1])$labels
  boot_med[b] <- adjustedRandIndex(lab, ref[idx])
}
cat(sprintf("bootstrap (median base): mean %.3f  95%% %.3f-%.3f\n",
            mean(boot_med), quantile(boot_med, 0.025), quantile(boot_med, 0.975)))

cat("\n--- stability summary (tableS3 is written by R/01_cohort_analyses.R) ---\n")
print(data.table(
  Comparison = c(
    sprintf("Median-imputation vs primary (mean ARI %.3f, range %.3f-%.3f across 100 re-inits)",
            mean(med_aris), min(med_aris), max(med_aris)),
    sprintf("Initialization stability, primary base (nstart=1, 25 re-inits, mean ARI %.3f, range %.3f-%.3f)",
            mean(init_aris), min(init_aris), max(init_aris)),
    sprintf("Bootstrap stability, primary base (500 x 70%% resamples, mean ARI %.3f, 95%% interval %.3f-%.3f)",
            mean(boot_mice), quantile(boot_mice, 0.025), quantile(boot_mice, 0.975)),
    sprintf("Bootstrap stability, median-imputation base (500 x 70%% resamples, mean ARI %.3f, 95%% interval %.3f-%.3f)",
            mean(boot_med), quantile(boot_med, 0.025), quantile(boot_med, 0.975))),
  ARI = c(mean(med_aris), mean(init_aris), mean(boot_mice), mean(boot_med)),
  n   = rep(nrow(Xm_s), 4L)))
cat("\n==== DONE ====\n")

# Threshold (15.3) stability in MIMIC subgroups
f <- fread(file.path(PROJECT_ROOT, "data", "feature_matrix_with_clusters.csv"))

run_stats <- function(y1, wbc, label) {
  ok <- !is.na(y1) & !is.na(wbc)
  y1 <- as.integer(y1[ok]); wbc <- wbc[ok]
  if (sum(y1) == 0 || sum(y1) == length(y1)) {
    return(data.table(Subgroup = label, n = length(y1), AUC = NA_real_,
                      Youden_cutoff = NA_real_, Sens_pct = NA_real_,
                      Spec_pct = NA_real_, PPV_pct = NA_real_, NPV_pct = NA_real_))
  }
  r <- roc(y1, wbc, quiet = TRUE)
  auc_v <- as.numeric(auc(r))
  yd <- r$sensitivities + r$specificities - 1
  cut <- r$thresholds[which.max(yd)]
  pred <- as.integer(wbc > 15.3)
  tp <- sum(pred == 1 & y1 == 1); fp <- sum(pred == 1 & y1 == 0)
  tn <- sum(pred == 0 & y1 == 0); fn <- sum(pred == 0 & y1 == 1)
  data.table(Subgroup = label, n = length(y1), AUC = auc_v,
             Youden_cutoff = cut,
             Sens_pct = tp / (tp + fn) * 100, Spec_pct = tn / (tn + fp) * 100,
             PPV_pct = tp / (tp + fp) * 100, NPV_pct = tn / (tn + fn) * 100)
}

y1 <- f$km_class == 1; wbc <- f$wbc_baseline
res <- run_stats(y1, wbc, "Overall (MIMIC)")
for (s in sort(unique(na.omit(f$infx_src)))) {
  idx <- f$infx_src == s
  res <- rbind(res, run_stats(y1[idx], wbc[idx], sprintf("Infection source: %s", s)))
}
res <- rbind(res, run_stats(y1[f$immunocompromised == 0],
                            wbc[f$immunocompromised == 0],
                            "Excluding immunocompromised"))
print(res)
fwrite(res, file.path(TABLES, "tableS7b_threshold_subgroups.csv"))
cat("\nSaved tableS7b_threshold_subgroups.csv\n")

# Finalize canonical tables (S6 / S10 / S12 / S13)

FEATS  <- c("wbc_baseline", "wbc_t24h", "wbc_t48h", "wbc_t72h")
CUTOFF <- 15.3

# Part A — baseline WBC distribution, AUC and outlier-exclusion robustness (S12)

df <- fread(file.path(PROJECT_ROOT, "data", "feature_matrix_with_clusters.csv"))
rows <- list()
add <- function(k, v) rows[[length(rows) + 1]] <<- data.table(Metric = k, Value = v)

# ---- A1. C1 baseline WBC distribution --------------------------------------
c1 <- df[km_class == 1, wbc_baseline]
c1 <- c1[!is.na(c1)]
nA  <- length(c1); mu <- mean(c1)
g1 <- mean((c1 - mu)^3) / (mean((c1 - mu)^2))^1.5      # biased sample skewness
skew_adj <- g1 * sqrt(nA * (nA - 1)) / (nA - 2)        # bias-corrected sample skewness

add("C1 n", sprintf("%d", nA))
add("C1 baseline WBC mean, K/uL", sprintf("%.1f", mu))
add("C1 baseline WBC median, K/uL", sprintf("%.1f", median(c1)))
add("C1 baseline WBC skewness", sprintf("%.2f", skew_adj))
add("C1 baseline WBC >40 n", sprintf("%d", sum(c1 > 40)))
add("C1 baseline WBC >40 pct", sprintf("%.1f%%", mean(c1 > 40) * 100))
add("C1 baseline WBC >50 n", sprintf("%d", sum(c1 > 50)))
add("C1 baseline WBC >50 pct", sprintf("%.1f%%", mean(c1 > 50) * 100))

# ---- A2. AUC with bootstrap CI + fixed-cutoff operating point --------------
ok <- !is.na(df$wbc_baseline)
yA <- as.integer(df$km_class[ok] == 1)
xA <- as.numeric(df$wbc_baseline[ok])

auc_rank <- function(yy, xx) {
  n1 <- sum(yy == 1); n0 <- sum(yy == 0)
  r <- rank(c(xx[yy == 1], xx[yy == 0]))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

auc_val <- auc_rank(yA, xA)
set.seed(20260910)
auc_bs <- vapply(seq_len(1000), function(i) {
  idx <- sample.int(length(yA), length(yA), replace = TRUE)
  if (length(unique(yA[idx])) < 2) return(NA_real_)
  auc_rank(yA[idx], xA[idx])
}, numeric(1))
auc_ci <- quantile(auc_bs, c(0.025, 0.975), na.rm = TRUE)

pred <- as.integer(xA > CUTOFF)
tp <- sum(pred == 1 & yA == 1); fp <- sum(pred == 1 & yA == 0)
tn <- sum(pred == 0 & yA == 0); fn <- sum(pred == 0 & yA == 1)

add("AUC baseline WBC -> C1", sprintf("%.3f", auc_val))
add("AUC 95% CI lower", sprintf("%.3f", auc_ci[[1]]))
add("AUC 95% CI upper", sprintf("%.3f", auc_ci[[2]]))
add("Threshold, K/uL", sprintf("%.1f", CUTOFF))
add("Sensitivity", sprintf("%.1f%%", tp / (tp + fn) * 100))
add("Specificity", sprintf("%.1f%%", tn / (tn + fp) * 100))
add("PPV", sprintf("%.1f%%", tp / (tp + fp) * 100))
add("NPV", sprintf("%.1f%%", tn / (tn + fn) * 100))

# ---- A3. Outlier-exclusion re-clustering -----------------------------------
X <- as.matrix(df[, ..FEATS]); mode(X) <- "numeric"
for (j in seq_len(ncol(X))) {
  na_idx <- is.na(X[, j]); if (any(na_idx)) X[na_idx, j] <- median(X[, j], na.rm = TRUE)
}
Xs <- scale(X); Xs[!is.finite(Xs)] <- 0

keepA <- df$wbc_baseline <= 40
set.seed(2026)
kmA <- kmeans(Xs[keepA, ], centers = 2, nstart = 50, iter.max = 100, algorithm = "Lloyd")
labA <- kmA$cluster
blA <- df$wbc_baseline[keepA]
if (mean(blA[labA == 1]) < mean(blA[labA == 2])) labA <- 3L - labA   # C1 = higher WBC
d28A <- df$death_28day[keepA]

add("Excl >40 re-cluster: n dropped", sprintf("%d", sum(!keepA)))
add("Excl >40 re-cluster: C1 mortality", sprintf("%.1f%%", mean(d28A[labA == 1]) * 100))
add("Excl >40 re-cluster: C2 mortality", sprintf("%.1f%%", mean(d28A[labA == 2]) * 100))

outA <- rbindlist(rows)
fwrite(outA, file.path(TABLES, "tableS12_baseline_distribution.csv"))
cat("-> tableS12_baseline_distribution.csv\n"); print(outA)

# Part B — K selection (S10) and three-cohort 15.3 threshold validation (S6)

# ---- B1. S10 = k_selection.csv + one K-selection note row ------------------
ksel <- file.path(TABLES, "k_selection.csv")
if (!file.exists(ksel)) stop("k_selection.csv missing - run R/01_cohort_analyses.R first")
k17 <- fread(ksel)
Ks <- k17$K[which.max(k17$Silhouette)]
Kg <- k17$K[which.max(ifelse(is.na(k17$Gap), -Inf, k17$Gap))]
k17_doc <- data.table(
  K = "—",
  WSS = "K selection criterion",
  Silhouette = sprintf("silhouette maximum at K = %g", Ks),
  Gap = sprintf("K = %g selected; gap statistic not used (its maximum lies at K = %g)", Ks, Kg))
s17 <- rbindlist(list(k17, k17_doc), fill = TRUE)
fwrite(s17, file.path(TABLES, "tableS10_k_selection.csv"))
cat("\n-> tableS10_k_selection.csv (k_selection.csv + one K-selection note row)\n")

# ---- B2. S6: 15.3 threshold performance across the three cohorts ----------
b1_path <- file.path(TABLES, "_b1_threshold_validation.csv")
if (!file.exists(b1_path))
  stop("_b1_threshold_validation.csv missing - run the small-analyses section first")
b1 <- fread(b1_path)

keepB <- b1$cohort == "MIMIC" |
         grepl("^eICU \\(de novo",  b1$cohort) |
         grepl("^NWICU \\(de novo", b1$cohort)
selB <- b1[keepB]

s6 <- data.table(
  Panel       = "A",
  Cohort      = c("MIMIC", "eICU", "NWICU"),
  Stratum     = "",
  n           = selB$n,
  AUC         = sprintf("%.4f", selB$auc),
  Cutoff      = sprintf("%.1f", selB$cutoff),
  Sensitivity = sprintf("%.1f", selB$sens_pct),
  Specificity = sprintf("%.1f", selB$spec_pct),
  PPV         = sprintf("%.1f", selB$ppv_pct),
  NPV         = sprintf("%.1f", selB$npv_pct),
  Mortality_pct = "", C1_pct = "")

# -- Panel B: mortality and C1 share by threshold --
b1m_path <- file.path(TABLES, "_b1_threshold_mortality.csv")
if (!file.exists(b1m_path))
  stop("_b1_threshold_mortality.csv missing - run the small-analyses section first")
b1m <- fread(b1m_path)
# wbc_hi == 1 is the WBC > 15.3 stratum
b1m[, Stratum := ifelse(wbc_hi == 1L, "WBC > 15.3", "WBC <= 15.3")]
b1m[, cohort := factor(cohort, levels = c("MIMIC", "eICU", "NWICU"))]
setorder(b1m, cohort, -wbc_hi)
s6 <- rbind(s6, data.table(
  Panel       = "B",
  Cohort      = as.character(b1m$cohort),
  Stratum     = b1m$Stratum,
  n           = b1m$n,
  AUC = "", Cutoff = "", Sensitivity = "", Specificity = "", PPV = "", NPV = "",
  Mortality_pct = sprintf("%.1f", b1m$mortality_pct),
  C1_pct        = sprintf("%.1f", b1m$c1_pct)), fill = TRUE)

.fm_wbc <- fread(file.path(PROJECT_ROOT, "data", "feature_matrix_with_clusters.csv"),
                 select = "wbc_baseline", showProgress = FALSE)$wbc_baseline
.mimic_eq <- sum(.fm_wbc == 15.3, na.rm = TRUE)
# Likewise the Youden optimum quoted in the note is read back from the table that
# computed it, instead of being written out as a literal.
.youden_raw <- fread(file.path(TABLES, "tableS7b_threshold_subgroups.csv"))[
  Subgroup == "Overall (MIMIC)"]$Youden_cutoff
s6 <- rbind(s6, data.table(
  Panel = "Note", Cohort = "", Stratum = "", n = "",
  AUC = "", Cutoff = "",
  Sensitivity = "Panel A metrics use WBC > 15.3 (not >=)",
  Specificity = sprintf("cases exactly at 15.3 (n = %d) fall on the low-risk side", .mimic_eq),
  PPV = sprintf("Youden optimum at %s", format(.youden_raw, trim = TRUE)),
  NPV = "Panel B reports mortality and C1 share within the same strata"), fill = TRUE)

fwrite(s6, file.path(TABLES, "tableS6_threshold_validation.csv"))
cat("-> tableS6_threshold_validation.csv\n"); print(s6)

# Part C — eICU time-anchor sensitivity (S13)

eicu_coh <- fread(file.path(PROJECT_ROOT, "data", "eicu_cohort.csv"))
eicu_asg <- fread(file.path(PROJECT_ROOT, "data", "eicu_cluster_assignment.csv"))

ea <- merge(eicu_asg[, .(patientunitstayid, km_class_recluster, hosp_death)],
            eicu_coh[, .(patientunitstayid, time_zero_min)],
            by = "patientunitstayid")
ea <- ea[!is.na(km_class_recluster) & !is.na(hosp_death)]
ea[, anchored := time_zero_min != 0]     # TRUE = genuine antibiotic time zero

grad <- function(d, lab) {
  m1 <- mean(d$hosp_death[d$km_class_recluster == 1]) * 100
  m2 <- mean(d$hosp_death[d$km_class_recluster == 2]) * 100
  data.table(
    Subset = lab,
    n = nrow(d),
    C1_n = sum(d$km_class_recluster == 1),
    C2_n = sum(d$km_class_recluster == 2),
    C1_mortality_pct = sprintf("%.1f%%", m1),
    C2_mortality_pct = sprintf("%.1f%%", m2),
    C1_C2_ratio = sprintf("%.2f", m1 / m2))
}

s20 <- rbindlist(list(
  grad(ea, "All eICU"),
  grad(ea[anchored == TRUE],  "Antibiotic-anchored (time zero != 0)"),
  grad(ea[anchored == FALSE], "ICU-admission fallback (time zero == 0)")))

fwrite(s20, file.path(TABLES, "tableS13_eicu_time_anchor_sensitivity.csv"))
cat("\n-> tableS13_eicu_time_anchor_sensitivity.csv\n"); print(s20)

cat("\n==== finalize tables DONE ====\n")
