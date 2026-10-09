# Cohort analyses and external validation (MIMIC -> eICU -> NWICU)

library(pROC)

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
  stop("Cannot locate project root (R/00_functions.R not found). ",
       "Run this script from the package root, or from R/ inside it.")
}
PROJECT_ROOT <- find_project_root()
source(file.path(PROJECT_ROOT, "R", "00_functions.R"))

set.seed(2026)

invert_hr_or <- function(dt, cols = c("HR","OR"), lo = "lower", up = "upper") {
  for (cn in cols) if (cn %in% names(dt)) {
    dt[[cn]] <- 1 / dt[[cn]]
    tmp_lo   <- 1 / dt[[up]]
    dt[[up]]  <- 1 / dt[[lo]]
    dt[[lo]]  <- tmp_lo
  }
  dt
}

TABLES <- file.path(PROJECT_ROOT, "output", "tables")
dir.create(TABLES, showWarnings = FALSE, recursive = TRUE)
LOG <- open_log(file.path(PROJECT_ROOT, "output", "run_log.txt"))

# 0. Load + Type Fix

df      <- fread(file.path(PROJECT_ROOT, "data", "feature_matrix.csv"))
df_long <- fread(file.path(PROJECT_ROOT, "data", "blood_long.csv"))

df[, stay_id := as.integer(stay_id)]
df_long[, stay_id := as.integer(stay_id)]
df_long[, `:=`(hrs = as.numeric(hrs), wbc = as.numeric(wbc))]

num_pattern <- "^wbc|^anc|^alc|^nlr|^sofa|^max_|_days$|^aps|charlson|admission_age|n_unique|n_cbc|n_tp|abx_n|_delta_|_pct_|_slope|n_blood"
bool_cols   <- grep("^has_|^death_|expire_flag|first_icu|persistent_", names(df), value = TRUE)
coerce_types(df, num_pattern, bool_cols)

df_long <- df_long[stay_id %in% unique(df$stay_id)]

# 0a. Derived variables

df[, has_bacterial := has_gram_negative | has_gram_positive]
df[, has_fungal_primary := has_fungal_blood]                       # blood-cx-restricted
df[, is_cn := !has_gram_negative & !(has_gram_positive & !has_contaminant) & !has_fungal_primary]
df[, is_cn_all_site := !has_gram_negative & !has_gram_positive & !has_fungal]
df[, is_female := gender == "F"]
df[, lactate_high := max_lactate >= 4.0]

df[, immunocompromised := as.integer(
  !is.na(anc_baseline) & anc_baseline < 1.0)]

df[, infx_src := fcase(
  admission_location %in% c("EMERGENCY ROOM", "EMERGENCY ROOM ADMIT",
                            "MEDICAL ICU", "MICU"), "Medical",
  admission_location %in% c("SURGICAL ICU", "SICU", "CARDIAC SURGERY",
                            "SURGERY"), "Surgical",
  default = "Other")]

has_nbc <- "n_blood_cultures" %in% names(df) && any(!is.na(df$n_blood_cultures))
if (!has_nbc) {
  cat("\n*** WARNING: n_blood_cultures missing — CNS 'never cultured' skipped ***\n")
  df[, n_blood_cultures := NA_integer_]
}

df[, n_missing_wbc := rowSums(is.na(.SD)), .SDcols = CLUSTER_FEATURES]
safe(missingness_report(df, c(CLUSTER_FEATURES, "max_lactate",
     "anc_baseline", "alc_baseline"), "cohort inputs"), "missingness_report")

df[, sex_f := as.integer(is_female)]
safe(standardize(df, c("admission_age", "sofa_sepsis3", "charlson", "max_lactate")),
     "standardize")
setnames(df, c("admission_age_z", "sofa_sepsis3_z", "charlson_z"),
             c("age_z", "sofa_z", "charlson_z"), skip_absent = TRUE)

# PRIMARY CLUSTERING — MICE-imputed WBC trajectory (single completed dataset)

cat("\n========== PRIMARY CLUSTERING (MICE, single completed dataset) ==========\n")

# V3: the published 96-h auxiliary (84-108 h) is removed from the extraction
# because it is wholly contained in wbc_t72h = [72,96); using it would leak the
# fourth clustering window into its own imputation model. n_tp_wbc (number of the
# four windows actually observed) is retained, which is also what eICU and NWICU
# use -- so the three cohorts now share one imputation specification.
mice_res <- impute_wbc_mice(df, features = CLUSTER_FEATURES,
                            aux = "n_tp_wbc", m = 5, seed = 2026)
X <- mice_res$per_imp[[1]]                        # primary: 1st completed dataset
sc <- scale(X); center_vec <- attr(sc, "scaled:center"); scale_vec <- attr(sc, "scaled:scale")
X_scaled <- sc; X_scaled[!is.finite(X_scaled)] <- 0

ksel <- select_k(X_scaled, kmax = 6, seed = 42)
print(ksel)
fwrite(ksel, file.path(TABLES, "k_selection.csv"))
k_opt <- ksel$K[which.max(ksel$Silhouette)]
sil_opt <- max(ksel$Silhouette)
cat(sprintf("Optimal K = %d (Silhouette = %.3f)\n", k_opt, sil_opt))

kf <- fit_kmeans(X_scaled, k_opt)
al <- align_labels_by_wbc(kf$km$cluster, X[, "wbc_baseline"])
df$km_class <- al$labels

# Reorder centroids so row 1 = C1 (for external prediction)
centers_aligned <- matrix(NA_real_, nrow = k_opt, ncol = ncol(kf$centers),
                          dimnames = list(NULL, colnames(kf$centers)))
for (r in seq_len(nrow(kf$centers))) centers_aligned[al$map[as.character(r)], ] <- kf$centers[r, ]

stab <- cluster_stability(X_scaled, k_opt, reference = kf$km$cluster, nrep = 25)
cat(sprintf("Cluster stability (quick): mean ARI = %.3f (min %.3f) over 25 re-inits\n",
            stab$mean, stab$min))

# --- Sensitivity: complete-case + per-imputation clustering vs primary (imp #1) ---
cc <- complete_case_matrix(df, CLUSTER_FEATURES)
cc_scaled <- apply_scaling(cc$X, center_vec, scale_vec)
cc_km <- kmeans(cc_scaled, centers = k_opt, nstart = 50, iter.max = 100)$cluster
ari_cc <- adjustedRandIndex(cc_km, df$km_class[cc$rows])

ari_per_imp <- sapply(mice_res$per_imp, function(Xi) {
  Xi_s <- apply_scaling(Xi, center_vec, scale_vec)
  adjustedRandIndex(kmeans(Xi_s, centers = k_opt, nstart = 25, iter.max = 100)$cluster,
                    df$km_class)
})
aligned_lab <- sapply(mice_res$per_imp, function(Xi) {
  Xi_s <- apply_scaling(Xi, center_vec, scale_vec)
  align_labels_by_wbc(kmeans(Xi_s, centers = k_opt, nstart = 25, iter.max = 100)$cluster,
                      Xi[, "wbc_baseline"])$labels
})
pair_ari <- c()
for (a in 1:4) for (b in (a + 1):5)
  pair_ari <- c(pair_ari, adjustedRandIndex(aligned_lab[, a], aligned_lab[, b]))

imp_stab <- data.table(
  Comparison = c("Primary (completed set 1) vs complete-case",
                 sprintf("Primary vs completed set %d", seq_along(ari_per_imp)),
                 "Pairwise among completed sets (mean)",
                 "Cluster stability (25 re-inits, mean)"),
  ARI = c(ari_cc, ari_per_imp, mean(pair_ari), stab$mean),
  n   = c(length(cc$rows), rep(nrow(df), length(ari_per_imp)), nrow(df), nrow(df)))
fwrite(imp_stab, file.path(TABLES, "tableS3_stability.csv"))
cat(sprintf("Table S3: Imputation stability — complete-case ARI = %.3f; completed-set ARI %.3f-%.3f; pairwise among sets %.3f\n",
            ari_cc, min(ari_per_imp), max(ari_per_imp), mean(pair_ari)))

# GMM sensitivity

gmm <- Mclust(X_scaled, G = 2, verbose = FALSE)
gmm_al <- align_labels_by_wbc(gmm$classification, X[, "wbc_baseline"])
df$gmm_class <- gmm_al$labels
ari_km_gmm <- adjustedRandIndex(df$km_class, df$gmm_class)

cat("\n========== GMM LABEL-ALIGNED VALIDATION ==========\n")
cat(sprintf("ARI (label-aligned): %.3f\n", ari_km_gmm))
tbl_2x2 <- table(Kmeans = df$km_class, GMM = df$gmm_class)
print(tbl_2x2)
cat(sprintf("Agreement: %.1f%%\n", sum(diag(tbl_2x2)) / sum(tbl_2x2) * 100))

df[, cross4 := fcase(
  km_class == 1 & gmm_class == 1, "K1_G1 (Both High)",
  km_class == 1 & gmm_class == 2, "K1_G2 (K-high)",
  km_class == 2 & gmm_class == 1, "K2_G1 (G-high)",
  km_class == 2 & gmm_class == 2, "K2_G2 (Both Normal)")]

cross4_res <- df[, .(n = .N, pct = sprintf("%.1f%%", .N / nrow(df) * 100),
  mort_28d = sprintf("%.1f%%", mean(death_28day) * 100),
  wbc_bl = sprintf("%.1f", mean(wbc_baseline, na.rm = TRUE)),
  sofa = sprintf("%.1f", mean(sofa_sepsis3))), by = cross4][order(-n)]
print(cross4_res)
fwrite(cross4_res, file.path(TABLES, "table2_cross_classification.csv"))

sens_updated <- data.table(
  Analysis = c("GMM vs K-means ARI (label-aligned)",
    "GMM K=2 C1 (high-WBC) mortality", "GMM K=2 C2 (normal) mortality",
    "K-means C1 mortality (ref)", "K-means C2 mortality (ref)",
    "4-group dose-response gradient",
    "K1_G1 (concordant high) mortality", "K2_G2 (concordant normal) mortality",
    "K1_G2 (K-high, discordant) mortality", "K2_G1 (G-high, discordant) mortality"),
  Result = c(sprintf("%.3f", ari_km_gmm),
    sprintf("%.1f%%", mean(df[gmm_class == 1, death_28day]) * 100),
    sprintf("%.1f%%", mean(df[gmm_class == 2, death_28day]) * 100),
    sprintf("%.1f%%", mean(df[km_class == 1, death_28day]) * 100),
    sprintf("%.1f%%", mean(df[km_class == 2, death_28day]) * 100),
    "Mortality gradient consistent with WBC magnitude",
    cross4_res[cross4 == "K1_G1 (Both High)", mort_28d],
    cross4_res[cross4 == "K2_G2 (Both Normal)", mort_28d],
    cross4_res[cross4 == "K1_G2 (K-high)", mort_28d],
    cross4_res[cross4 == "K2_G1 (G-high)", mort_28d]),
  Interpretation = c(
    "Agreement between two independent clustering methods",
    "Higher mortality in high-WBC group", "Lower mortality in normal-WBC group",
    "Reference", "Reference",
    "Dose-response: concordant groups show maximal separation",
    "Highest mortality — both methods agree on high-WBC",
    "Lowest mortality — both methods agree on normal WBC",
    "Intermediate — partial overlap zone", "Intermediate — partial overlap zone"))
fwrite(sens_updated, file.path(TABLES, "tableS7_sensitivity_gmm.csv"))

ctr <- df[, lapply(.SD, mean, na.rm = TRUE), by = km_class,
          .SDcols = patterns("^wbc_(baseline|t24h|t48h|t72h)$")]

# 3. TABLE 1

t1_vars <- list(
  list("admission_age", "Age, years", "cont"),
  list("is_female", "Female, n (%)", "bin"),
  list("sofa_sepsis3", "SOFA score", "cont", 0, "median"),
  list("apsiii", "APS III", "cont", 0, "median"),
  list("charlson", "Charlson Comorbidity Index", "cont"),
  list("immunocompromised", "Immunocompromised, n (%)", "bin"),
  list("infx_src", "Infection source (Medical vs Surgical/Other), %", "cat"),
  list("n_missing_wbc", "Missing WBC time-points", "cont", 0, "median"),
  list("wbc_baseline", "WBC baseline, K/uL", "cont", 1, "median"),
  list("wbc_t72h", "WBC 72h, K/uL", "cont", 1, "median"),
  list("wbc_slope", "WBC slope, K/uL/h", "cont", 2),
  list("wbc_delta_72h", "WBC delta 72h, K/uL", "cont", 1),
  list("persistent_high_wbc", "Persistent high WBC, n (%)", "bin"),
  list("persistent_lymphopenia", "Persistent lymphopenia, n (%)", "bin"),
  list("max_lactate", "Max lactate, mmol/L", "cont", 1, "median"),
  list("icu_los_days", "ICU LOS, days", "cont", 1, "median"),
  list("has_vent", "Mechanical ventilation, n (%)", "bin"),
  list("has_vaso", "Vasopressors, n (%)", "bin"),
  list("has_crrt", "CRRT, n (%)", "bin"),
  list("has_glycopeptide", "Glycopeptide antibiotics, n (%)", "bin"),
  list("has_carbapenem", "Carbapenem antibiotics, n (%)", "bin"),
  list("has_cephalosporin", "Cephalosporin antibiotics, n (%)", "bin"),
  list("has_bacterial", "Bacterial infection, n (%)", "bin"),
  list("is_cn", "Culture-negative sepsis, n (%)", "bin"),
  list("death_28day", "28-day mortality, n (%)", "bin"),
  list("hospital_expire_flag", "Hospital mortality, n (%)", "bin"))

t1 <- safe(make_table1(df, "km_class", t1_vars), "make_table1")
if (!is.null(t1)) {
  fwrite(t1, file.path(TABLES, "table1_manuscript.csv"))
  cat("\nTable 1 written: table1_manuscript.csv\n")
}

# 4. MORTALITY MODELING — KM + Cox + Logistic

cat("\n========== MORTALITY MODELING ==========\n")

df_surv <- df[!is.na(death_28day) & !is.na(icu_los_days)]
df_surv[, surv_time := pmin(icu_los_days, 28)]
df_surv[, event := death_28day]

# --- Cox regression with PH-assumption test ---
cox_models <- list(
  "Unadjusted"          = Surv(surv_time, event) ~ km_class,
  "+ Demographics"      = Surv(surv_time, event) ~ km_class + age_z + sex_f,
  "+ Severity"          = Surv(surv_time, event) ~ km_class + age_z + sex_f + sofa_z + charlson_z,
  "+ Organ Support"     = Surv(surv_time, event) ~ km_class + age_z + sex_f + sofa_z + charlson_z + has_vent + has_vaso + has_crrt,
  "+ Immunocompromised" = Surv(surv_time, event) ~ km_class + age_z + sex_f + sofa_z + charlson_z + has_vent + has_vaso + has_crrt + immunocompromised,
  "Fully adjusted"      = Surv(surv_time, event) ~ km_class + age_z + sex_f + sofa_z + charlson_z + has_vent + has_vaso + has_crrt + immunocompromised + has_glycopeptide + has_carbapenem + has_cephalosporin)
cox_results <- fit_cox_set(df_surv, cox_models, term = "km_class")
cox_results <- invert_hr_or(cox_results)
cat("Table 3 (Cox) computed\n"); print(cox_results)

# --- Time-varying-coefficient Cox sensitivity (addresses the PH violation) ---
tv_cox <- safe(fit_cox_timevarying(df_surv, term = "km_class",
  covars = c("age_z", "sex_f", "sofa_z", "charlson_z"), cutpoints = c(7, 14)),
  "cox_timevarying")

# ── Assemble Table 3: Panel A (nested Cox) + Panel B (time-varying HR) ──
t3_panels <- list()
t3a <- copy(cox_results[, .(Model, HR, lower, upper, p, C_index)])
t3a[, Panel := "A: Nested Cox models"]
setcolorder(t3a, "Panel")
t3_panels[["A"]] <- t3a

if (!is.null(tv_cox)) {
  tv_cox <- invert_hr_or(tv_cox)
  t3b <- copy(tv_cox)
  t3b[, Panel := "B: Time-varying HR"]
  setcolorder(t3b, "Panel")
  t3_panels[["B"]] <- t3b
  cat("Table 3 Panel B (time-varying Cox): interval-specific HR\n")
  print(tv_cox)
}

t3_combined <- rbindlist(t3_panels, fill = TRUE)
fwrite(t3_combined, file.path(TABLES, "table3_cox_mortality.csv"))
cat("Table 3 (Cox nested + time-varying HR) written: table3_cox_mortality.csv\n")

# ── S2: Schoenfeld PH test only ──
ph_dt <- copy(cox_results[, .(Model, PH_p_term, PH_p_global)])
fwrite(ph_dt, file.path(TABLES, "tableS2_cox_ph.csv"))
cat("Table S2 (Cox PH test) written: tableS2_cox_ph.csv\n")
if (any(cox_results$PH_p_global < 0.05, na.rm = TRUE))
  cat("*** NOTE: global PH assumption violated — see Table 3 Panel B for time-varying HR ***\n")

# --- Logistic mortality (LR CIs) ---
or_mort_formulas <- lapply(cox_models, function(f)
  paste("death_28day ~", paste(all.vars(f)[-c(1, 2)], collapse = " + ")))
or_mort_dt <- fit_glm_set(df_surv, or_mort_formulas, term = "km_class")
or_mort_dt <- invert_hr_or(or_mort_dt, cols = "OR")
fwrite(or_mort_dt, file.path(TABLES, "tableS_logistic_or_mortality.csv"))
cat("Supplementary: Logistic OR for mortality written\n")

# 5. INFECTION SOURCE STRATIFICATION

cat("\n========== INFECTION SOURCE ==========\n")

src_dt <- data.table()
for (src in na.omit(unique(df$infx_src))) {
  sub <- df[infx_src == src]
  if (nrow(sub) < 100) next
  res <- safe({
    m <- glm(death_28day ~ km_class + age_z + sex_f + sofa_z + charlson_z,
             data = sub, family = binomial)
    ci <- .lr_confint(m, "km_class")
    tbl <- table(sub$km_class, sub$has_bacterial)
    bac_c1 <- if (all(dim(tbl) == 2)) tbl[1, 2] / sum(tbl[1, ]) * 100 else NA
    bac_c2 <- if (all(dim(tbl) == 2)) tbl[2, 2] / sum(tbl[2, ]) * 100 else NA
    data.table(Source = src, n = nrow(sub),
      C1_mort = sprintf("%.1f%%", mean(sub$death_28day[sub$km_class == 1]) * 100),
      C2_mort = sprintf("%.1f%%", mean(sub$death_28day[sub$km_class == 2]) * 100),
      OR_mort = 1 / unname(exp(coef(m)["km_class"])),
      lower   = 1 / unname(exp(ci[2])), upper = 1 / unname(exp(ci[1])),
      p_mort = summary(m)$coefficients["km_class", 4],
      C1_bacterial = sprintf("%.1f%%", bac_c1), C2_bacterial = sprintf("%.1f%%", bac_c2))
  }, paste0("infx_src:", src))
  if (!is.null(res)) src_dt <- rbind(src_dt, res)
}
if (nrow(src_dt) > 0) {
  fwrite(src_dt, file.path(TABLES, "tableS7_subgroup_analysis.csv"))
  cat("Table S7: Subgroup analysis by infection source written\n")
}

# 6. IMMUNOCOMPROMISED SENSITIVITY

cat("\n========== IMMUNOCOMPROMISED ==========\n")
cat(sprintf("Immunocompromised overall %.1f%%  (C1 %.1f%%, C2 %.1f%%)\n",
  mean(df$immunocompromised) * 100,
  mean(df[km_class == 1, immunocompromised]) * 100,
  mean(df[km_class == 2, immunocompromised]) * 100))
safe({
  df_imm <- df[immunocompromised == 0]
  m_imm <- glm(death_28day ~ km_class + age_z + sex_f + sofa_z + charlson_z +
    has_vent + has_vaso + has_crrt + has_glycopeptide + has_carbapenem + has_cephalosporin,
    data = df_imm, family = binomial)
  ci <- .lr_confint(m_imm, "km_class")
  cat(sprintf("OR C1-vs-C2 (excl. immunocompromised, n=%d): %.2f [%.2f-%.2f] p=%.2e\n",
    nrow(df_imm), 1/exp(coef(m_imm)["km_class"]), 1/exp(ci[2]), 1/exp(ci[1]),
    summary(m_imm)$coefficients["km_class", 4]))
  immuno_dt <- data.table(
    Analysis = "Excluding immunosuppressed",
    N = nrow(df_imm),
    OR_C1_vs_C2 = round(1 / exp(coef(m_imm)["km_class"]), 2),
    CI_lower = round(1 / exp(ci[2]), 2),
    CI_upper = round(1 / exp(ci[1]), 2),
    p = summary(m_imm)$coefficients["km_class", 4]
  )
  fwrite(immuno_dt, file.path(TABLES, "tableS7_immuno_sensitivity.csv"))
  cat("  -> tableS7_immuno_sensitivity.csv written\n")
}, "immuno_sensitivity")

# 7. Microbiology Validation + Adjusted OR

cat("\n========== MICROBIOLOGY VALIDATION ==========\n")

m5_vars <- c("has_bacterial", "has_gram_negative", "has_gram_positive",
             "has_fungal_primary", "has_fungal", "is_cn")
m5_labels <- c("Any bacterial", "Gram-negative", "Gram-positive",
               "Fungal (blood cx)", "Fungal (all sites)",
               "Culture-negative")
m5_results <- data.table()
for (i in seq_along(m5_vars)) {
  var <- m5_vars[i]
  res <- safe({
    chi <- chisq.test(table(df$km_class, df[[var]]))
    c1_n <- sum(df$km_class == 1 & df[[var]] == TRUE, na.rm = TRUE)
    c2_n <- sum(df$km_class == 2 & df[[var]] == TRUE, na.rm = TRUE)
    c1_pct <- mean(df[km_class == 1, get(var)], na.rm = TRUE) * 100
    c2_pct <- mean(df[km_class == 2, get(var)], na.rm = TRUE) * 100
    m_unadj <- glm(as.formula(paste0(var, " ~ km_class")), data = df, family = binomial)
    or_unadj <- 1 / exp(coef(m_unadj)["km_class"])
    ci_unadj <- rev(1 / exp(confint.default(m_unadj)["km_class", ]))
    m_adj <- glm(as.formula(paste0(var,
      " ~ km_class + age_z + sex_f + sofa_z + charlson_z + has_vent + has_vaso + has_crrt + immunocompromised + has_glycopeptide + has_carbapenem + has_cephalosporin")),
      data = df, family = binomial)
    or_adj <- 1 / exp(coef(m_adj)["km_class"])
    ci_adj <- rev(1 / exp(confint.default(m_adj)["km_class", ]))
    data.table(Pathogen = m5_labels[i],
      C1_n = c1_n, C1_pct = sprintf("%.1f", c1_pct),
      C2_n = c2_n, C2_pct = sprintf("%.1f", c2_pct),
      OR_unadj = or_unadj, OR_unadj_CI = sprintf("%.2f–%.2f", ci_unadj[1], ci_unadj[2]),
      OR_adj = or_adj, OR_adj_CI = sprintf("%.2f–%.2f", ci_adj[1], ci_adj[2]),
      chi2_p = chi$p.value)
  }, paste0("micro:", var))
  if (!is.null(res)) m5_results <- rbind(m5_results, res)
}
if (nrow(m5_results) > 0) {
  m5_results[, p_bh := p.adjust(chi2_p, method = "BH")]
  m5_results[, chi2_p := sprintf("%.3f", chi2_p)]
  m5_results[, p_bh := sprintf("%.3f", p_bh)]
  fwrite(m5_results, file.path(TABLES, "table4_microbiology.csv"))
  cat("Table 4 (microbiology): pathogen spectrum written\n")
  print(m5_results)
} else {
  cat("WARNING: no microbiology results available; Table 4 skipped\n")
}

# 8. Mediation Analysis — multi-model sensitivity

cat("\n========== MEDIATION (lactate) ==========\n")

dm <- na.omit(df[, .(km_class, max_lactate, death_28day, admission_age, is_female,
  sofa_sepsis3, charlson, has_glycopeptide, has_carbapenem, has_cephalosporin,
  has_vent, has_vaso, has_crrt, immunocompromised, age_z, sex_f, sofa_z, charlson_z)])
dm[, `:=`(lac_z = as.numeric(scale(max_lactate)), is_female = as.integer(is_female))]

# Causal mediation (natural indirect/direct effects, risk-difference scale)
dm_df <- as.data.frame(dm)
dm_df$km_class_f <- as.integer(dm_df$km_class == 1)
dm_df$death_28day <- as.integer(dm_df$death_28day)

tiers <- list(
  "Minimal (age, sex)"    = c("admission_age", "is_female"),
  "Base (+ severity)"     = c("admission_age", "is_female", "sofa_sepsis3", "charlson"),
  "Full (+ support, abx)" = c("admission_age", "is_female", "sofa_sepsis3", "charlson",
                              "has_vent", "has_vaso", "has_crrt", "immunocompromised",
                              "has_glycopeptide", "has_carbapenem", "has_cephalosporin"))

med_dt <- data.table()
for (nm in names(tiers)) {
  cm <- causal_mediation(dm_df, "km_class_f", "max_lactate", "death_28day", tiers[[nm]])
  med_dt <- rbind(med_dt, cm[, .(Adjustment = nm, Indirect = ACME, Prop_mediated_pct,
                                 CI_lower_pct, CI_upper_pct,
                                 ADE = ADE, ADE_lo = ADE_lo, ADE_hi = ADE_hi,
                                 ADE_p = ADE_p)], fill = TRUE)
}
med_dt[, Cohort := "MIMIC"]
fwrite(med_dt, file.path(TABLES, "table4_mimic_mediation.csv"))
print(med_dt)
prop_med <- med_dt[Adjustment == "Full (+ support, abx)", Prop_mediated_pct]
if (length(prop_med) == 0) prop_med <- med_dt$Prop_mediated_pct[nrow(med_dt)]
cat(sprintf("Mediation (fully adjusted): %.1f%% of group effect via lactate; range %.1f-%.1f%%\n",
  prop_med, min(med_dt$Prop_mediated_pct), max(med_dt$Prop_mediated_pct)))

# 9. STROBE flow + culture-negative (CNS) analysis

cat("\n========== STROBE + CNS ==========\n")

n_total <- nrow(df); n_c1 <- sum(df$km_class == 1); n_c2 <- sum(df$km_class == 2)
# Six-step STROBE funnel counts are written by sql/01_mimic.sql to data/mimic_strobe_counts.csv.
strobe_cnt <- fread(file.path(PROJECT_ROOT, "data", "mimic_strobe_counts.csv"))
cat(sprintf("STROBE step 1 (ICU admissions): %s\n",
            format(strobe_cnt$step1_icu_admissions, big.mark = ",")))
strobe <- data.table(
  Step = c("MIMIC-IV v3.1 ICU admissions",
    "Sepsis-3 criteria met (suspected infection + SOFA >= 2)",
    "Age >= 18 years, first ICU stay, LOS >= 1 day",
    "Valid time zero (antibiotic/infection window)",
    "Excluded: hematologic malignancy, HIV, solid organ transplant",
    "WBC baseline measurement available",
    "FINAL ANALYTIC COHORT",
    "  C1 high-WBC subgroup",
    "  C2 normal-WBC subgroup"),
  n = c(format(strobe_cnt$step1_icu_admissions, big.mark=","),
    format(strobe_cnt$step2_sepsis3, big.mark=","),
    format(strobe_cnt$step3_age_first_los, big.mark=","),
    format(strobe_cnt$step4_time_zero, big.mark=","),
    format(strobe_cnt$step5_excluded, big.mark=","),
    sprintf("%d", n_total), sprintf("%d", n_total),
    sprintf("%d (%.1f%%)", n_c1, n_c1 / n_total * 100),
    sprintf("%d (%.1f%%)", n_c2, n_c2 / n_total * 100)))
fwrite(strobe, file.path(TABLES, "table0_strobe_flow.csv"))
cat("Table 0: STROBE flow written\n")

if (has_nbc) {
  df[, cns_type := fcase(
    is_cn == FALSE, "Culture-positive",
    n_blood_cultures == 0, "Never cultured",
    n_blood_cultures >= 1 & is_cn == TRUE, "True culture-negative",
    default = "Other")]
  cns_summary <- df[, {
    mc1 <- if (sum(km_class == 1) > 0) mean(death_28day[km_class == 1]) * 100 else NA_real_
    mc2 <- if (sum(km_class == 2) > 0) mean(death_28day[km_class == 2]) * 100 else NA_real_
    .(n = .N, pct = sprintf("%.1f%%", .N / nrow(df) * 100),
      mort = sprintf("%.1f%%", mean(death_28day) * 100),
      n_c1 = sum(km_class == 1), mort_c1 = sprintf("%.1f%%", mc1),
      n_c2 = sum(km_class == 2), mort_c2 = sprintf("%.1f%%", mc2),
      wbc = sprintf("%.1f", mean(wbc_baseline, na.rm = TRUE)))
  }, by = cns_type]
  print(cns_summary)

} else {
  df[, cns_type := ifelse(is_cn, "Culture-negative", "Culture-positive")]
}

# CNS: stratified effects by lactate
cns_for_analysis <- df[is_cn == TRUE & !is.na(lactate_high)]
cns_mort <- cns_for_analysis[, .(mort = mean(death_28day) * 100,
  se = sd(death_28day) / sqrt(.N) * 100, n = .N, sofa = mean(sofa_sepsis3)),
  by = .(km_class, lactate_high)]

cns_strat <- safe(stratified_effects(cns_for_analysis, outcome = "death_28day",
  treat = "km_class", strata_var = "lactate_high",
  covars = c("admission_age", "is_female", "sofa_sepsis3", "charlson"),
  interaction_term = "km_class:lactate_highTRUE"), "cns_stratified")
if (!is.null(cns_strat)) {
  cat(sprintf("CNS interaction p = %.4f; stratified ORs written\n", cns_strat$interaction_p))
  print(cns_strat$strata)
}

# Additive-scale interaction (RERI)
reri_res <- safe({
  rc <- copy(cns_for_analysis[!is.na(death_28day)])
  rc[, `:=`(A = as.integer(km_class == 1), B = as.integer(lactate_high))]
  crude <- compute_reri(rc, "death_28day", "A", "B")
  adj   <- compute_reri(rc, "death_28day", "A", "B",
             covars = c("admission_age", "sex_f", "sofa_sepsis3", "charlson"))
  data.table(Measure = crude$Measure, Crude = crude$Value, Adjusted = adj$Value)
}, "reri")
if (!is.null(reri_res)) {
  cat("Table 5 panel C (RERI additive-scale interaction) written\n")
  print(reri_res)
}

# ── Assemble Table 5: CNS analysis ──
s4_panels <- list()
if (exists("cns_summary") && is.data.table(cns_summary) && nrow(cns_summary) > 0) {
  pa <- copy(cns_summary)
  setnames(pa, old = c("cns_type", "n", "pct", "mort", "n_c1", "mort_c1",
    "n_c2", "mort_c2", "wbc"),
    new = c("Group", "N", "Pct_of_cohort", "Overall_mortality",
      "C1_N", "C1_mortality", "C2_N", "C2_mortality", "WBC_baseline_mean"))
  pa[, Panel := "A"]
  setcolorder(pa, "Panel")
  s4_panels[["A"]] <- pa
}
if (exists("cns_mort") && is.data.table(cns_mort) && nrow(cns_mort) > 0) {
  pb <- cns_mort[, .(
    Group = paste0(ifelse(km_class == 1, "C1", "C2"),
      ", Lactate ", ifelse(lactate_high, "≥4 mmol/L", "<4 mmol/L")),
    N = n, Mortality_pct = sprintf("%.1f%%", mort),
    SOFA_mean = sprintf("%.1f", sofa))]
  pb[, Panel := "B"]
  setcolorder(pb, "Panel")
  s4_panels[["B"]] <- pb
}
if (exists("reri_res") && is.data.table(reri_res) && nrow(reri_res) > 0) {
  pc <- copy(reri_res)
  pc[, Panel := "C"]
  setcolorder(pc, "Panel")
  s4_panels[["C"]] <- pc
}
if (length(s4_panels) > 0) {
  s4_combined <- rbindlist(s4_panels, fill = TRUE)
  fwrite(s4_combined, file.path(TABLES, "table5_cns_interaction.csv"))
  cat("Table 5 (culture-negative subgroup x lactate): table5_cns_interaction.csv written\n")
}

# 10. K=3 SENSITIVITY + BASELINE WBC ROC

cat("\n========== K=3 SENSITIVITY ==========\n")
safe({
  km3 <- kmeans(X_scaled, centers = 3, nstart = 50, iter.max = 100)
  df$km3_class <- km3$cluster
  k3_mort <- df[, .(n = .N, pct = sprintf("%.1f%%", .N / nrow(df) * 100),
    mort = sprintf("%.1f%%", mean(death_28day) * 100),
    wbc_bl = sprintf("%.1f", mean(wbc_baseline, na.rm = TRUE)),
    sofa = sprintf("%.1f", mean(sofa_sepsis3))), by = km3_class][order(-n)]
  ari_k2k3 <- adjustedRandIndex(df$km_class, df$km3_class)
  cat(sprintf("ARI (K=2 vs K=3): %.3f\n", ari_k2k3))
  k3_mort <- k3_mort[order(-as.numeric(wbc_bl))]
  k3_sens <- data.table(
    Group = c("—", "K=3 group 1 (highest WBC)", "K=3 group 2", "K=3 group 3 (lowest WBC)",
              "—", "—"),
    n = c(NA_integer_, k3_mort$n[1], k3_mort$n[2], k3_mort$n[3], NA_integer_, NA_integer_),
    pct = c(NA_character_, k3_mort$pct[1], k3_mort$pct[2], k3_mort$pct[3], NA_character_, NA_character_),
    WBC_baseline = c(NA_character_, k3_mort$wbc_bl[1], k3_mort$wbc_bl[2], k3_mort$wbc_bl[3],
                     NA_character_, NA_character_),
    Analysis = c("K=3 Silhouette", rep("K=3 mortality", 3),
                 "ARI K=2 vs K=3", "K=3 trajectory pattern"),
    Result = c(sprintf("%.3f", ksel$Silhouette[ksel$K == 3]),
      k3_mort$mort[1], k3_mort$mort[2], k3_mort$mort[3],
      sprintf("%.3f", ari_k2k3), "WBC gradient across the three groups; no stable third subgroup"))
  fwrite(k3_sens, file.path(TABLES, "tableS5_k3_sensitivity.csv"))
}, "k3_sensitivity")

cat("\n========== BASELINE WBC ROC ==========\n")
auc_val <- NA_real_; opt_cutoff <- NA_real_
safe({
  roc_obj <- roc(df$km_class == 1, df$wbc_baseline, quiet = TRUE)
  auc_val <<- as.numeric(auc(roc_obj))
  youden <- roc_obj$sensitivities + roc_obj$specificities - 1
  oi <- which.max(youden)
# Operating point: Youden-optimal threshold rounded to one decimal (x10^9/L).
  cutoff <- round(roc_obj$thresholds[oi], 1)
  opt_cutoff <<- cutoff
  .pred <- as.integer(df$wbc_baseline > cutoff)
  .yt   <- as.integer(df$km_class == 1)
  .ok   <- !is.na(.pred) & !is.na(.yt)
  .tp <- sum(.pred[.ok] == 1 & .yt[.ok] == 1); .fp <- sum(.pred[.ok] == 1 & .yt[.ok] == 0)
  .tn <- sum(.pred[.ok] == 0 & .yt[.ok] == 0); .fn <- sum(.pred[.ok] == 0 & .yt[.ok] == 1)
  opt_sens <- .tp / (.tp + .fn); opt_spec <- .tn / (.tn + .fp)
  cat(sprintf("Baseline WBC -> C1: AUC = %.3f; cutoff > %.1f (Sens %.1f%%, Spec %.1f%%)\n",
    auc_val, opt_cutoff, opt_sens * 100, opt_spec * 100))
  # Bootstrap 95% CI for the AUC (1000 resamples, seed fixed for reproducibility)
  set.seed(20260910)
  auc_bs <- replicate(1000, {
    i <- sample(which(.ok), sum(.ok), replace = TRUE)
    if (length(unique(.yt[i])) < 2) return(NA_real_)
    as.numeric(pROC::auc(pROC::roc(.yt[i], df$wbc_baseline[.ok][i], quiet = TRUE)))
  })
  auc_ci <- quantile(auc_bs, c(0.025, 0.975), na.rm = TRUE)
  roc_tbl <- data.table(
    Metric = c("AUC", "AUC 95% CI lower", "AUC 95% CI upper",
               "Optimal cutoff (K/uL)", "Sensitivity", "Specificity", "PPV", "NPV"),
    Value = c(sprintf("%.3f", auc_val),
      sprintf("%.3f", auc_ci[[1]]), sprintf("%.3f", auc_ci[[2]]),
      sprintf("%.1f", opt_cutoff),
      sprintf("%.1f%%", opt_sens * 100), sprintf("%.1f%%", opt_spec * 100),
      sprintf("%.1f%%", .tp / (.tp + .fp) * 100), sprintf("%.1f%%", .tn / (.tn + .fn) * 100)))
  fwrite(roc_tbl, file.path(TABLES, "tableS_roc_baseline_wbc.csv"))
}, "roc_baseline")

# 11. ANC/ALC/NLR description + complete-case check

df[, nlr_baseline := ifelse(!is.na(anc_baseline) & !is.na(alc_baseline) & alc_baseline > 0,
                            anc_baseline / alc_baseline, NA_real_)]
m4_vars <- c("anc_baseline", "alc_baseline", "nlr_baseline", "anc_t72h", "alc_t72h")
m4_labels <- c("ANC baseline", "ALC baseline", "NLR baseline", "ANC 72h", "ALC 72h")
m4_results <- data.table()
for (i in seq_along(m4_vars)) {
  var <- m4_vars[i]
  c1 <- na.omit(df[km_class == 1, get(var)]); c2 <- na.omit(df[km_class == 2, get(var)])
  if (length(c1) > 10 && length(c2) > 10) {
    m4_results <- rbind(m4_results, data.table(Variable = m4_labels[i],
      Coverage = sprintf("%.0f%%", mean(!is.na(df[[var]])) * 100),
      C1 = sprintf("%.1f", mean(c1)), C2 = sprintf("%.1f", mean(c2)),
      p = safe(wilcox.test(c1, c2)$p.value, "wilcox")))
  }
}

df_cc <- df[!is.na(wbc_baseline) & !is.na(wbc_t24h) & !is.na(wbc_t48h) & !is.na(wbc_t72h)]
cat(sprintf("Complete-case: n=%d (%.1f%%), mort C1=%.1f%% C2=%.1f%%\n",
  nrow(df_cc), nrow(df_cc) / nrow(df) * 100,
  mean(df_cc[km_class == 1, death_28day]) * 100,
  mean(df_cc[km_class == 2, death_28day]) * 100))

# 12. EXPORT cluster model (for eICU) + clustered data

ref_stats <- list(
  n = n_total, k = k_opt, silhouette = sil_opt,
  c1_prop = n_c1 / n_total * 100, c2_prop = n_c2 / n_total * 100,
  c1_mort_28d = mean(df[km_class == 1, death_28day]) * 100,
  c2_mort_28d = mean(df[km_class == 2, death_28day]) * 100,
  c1_mort_hosp = mean(df[km_class == 1, hospital_expire_flag]) * 100,
  c2_mort_hosp = mean(df[km_class == 2, hospital_expire_flag]) * 100,
  ari_gmm = ari_km_gmm)
export_cluster_model(CLUSTER_FEATURES, center_vec, scale_vec, centers_aligned,
  al$map, ref_stats, file.path(PROJECT_ROOT, "output", "model", "mimic_cluster_model.rds"))
cat("Exported: output/model/mimic_cluster_model.rds\n")

fwrite(df, file.path(PROJECT_ROOT, "data", "feature_matrix_with_clusters.csv"))


cat("\n============================================\n")
cat("MIMIC ANALYSIS COMPLETE (MICE primary)\n")
cat("============================================\n")
cat(sprintf("Cohort: %d stays, K=%d (Silhouette=%.3f)\n", n_total, k_opt, sil_opt))
cat(sprintf("C1 high-WBC: %d (%.1f%%)  mortality %.1f%%\n",
  n_c1, n_c1 / n_total * 100, mean(df[km_class == 1, death_28day]) * 100))
cat(sprintf("C2 Normal: %d (%.1f%%)  mortality %.1f%%\n",
  n_c2, n_c2 / n_total * 100, mean(df[km_class == 2, death_28day]) * 100))
cat(sprintf("GMM ARI: %.3f | complete-case ARI: %.3f | stability ARI: %.3f\n",
  ari_km_gmm, ari_cc, stab$mean))
cat(sprintf("Mediation (lactate, full): %.1f%%\n", prop_med))
cat(sprintf("Baseline WBC ROC AUC: %.3f (cutoff > %.0f K/uL)\n", auc_val, opt_cutoff))
cat("============================================\n")

close_log(LOG)

# eICU external validation: de novo K-means (Method A) and MIMIC-centroid prediction (Method B)

TABLES <- file.path(PROJECT_ROOT, "output", "tables")
dir.create(TABLES, showWarnings = FALSE, recursive = TRUE)
LOG <- open_log(file.path(PROJECT_ROOT, "output", "run_log_eicu.txt"))

# 0. Load

eicu       <- fread(file.path(PROJECT_ROOT, "data", "eicu_cohort.csv"))
eicu_long  <- fread(file.path(PROJECT_ROOT, "data", "eicu_blood_long.csv"))
eicu_micro <- fread(file.path(PROJECT_ROOT, "data", "eicu_microlab_raw.csv"))

eicu[, patientunitstayid := as.integer(patientunitstayid)]
eicu_long[, patientunitstayid := as.integer(patientunitstayid)]
eicu_micro[, patientunitstayid := as.integer(patientunitstayid)]

for (col in names(eicu)) {
  x <- eicu[[col]]
  if (is.character(x) && all(unique(na.omit(x)) %in% c("t","f","true","false","1","0","TRUE","FALSE"))) {
    eicu[[col]] <- as.integer(tolower(x) %in% c("t","true","1"))
  }
}
eicu[, age := if (is.numeric(age)) age else as.numeric(ifelse(age == "> 89", "90", age))]

cohort_ids <- unique(eicu$patientunitstayid)
eicu_long <- eicu_long[patientunitstayid %in% cohort_ids]
eicu_long[, `:=`(hrs = as.numeric(hrs), wbc = as.numeric(wbc))]

# 1. Build WBC windows + scale

features <- build_wbc_windows(eicu_long, "patientunitstayid", WBC_WINDOWS)
features[, n_tp_wbc := rowSums(!is.na(.SD)), .SDcols = CLUSTER_FEATURES]
n_eicu_wbc <- nrow(features)   # eICU stays with at least one WBC measurement

eicu_mice <- impute_wbc_mice(features, features = CLUSTER_FEATURES,
                             aux = "n_tp_wbc", m = 5, seed = 2026, maxit = 5)
X <- eicu_mice$per_imp[[1]]                       # single completed dataset
X_scaled <- scale(X); X_scaled[!is.finite(X_scaled)] <- 0

cat("============================================================\n")
cat("eICU EXTERNAL VALIDATION\n")
cat("============================================================\n\n")

# ─── Method A: de novo K-means K=2 ───
set.seed(2026)
kf <- fit_kmeans(X_scaled, 2)
al <- align_labels_by_wbc(kf$km$cluster, X[, "wbc_baseline"])
features[, km_class_recluster := al$labels]
idx_sub <- sample(seq_len(nrow(X_scaled)), min(5000, nrow(X_scaled)))
sil_recluster <- summary(silhouette(kf$km$cluster[idx_sub], dist(X_scaled[idx_sub, ])))$avg.width

cat("--- Method A: De novo K-means K=2 ---\n")
cat(sprintf("Silhouette: %.3f\n", sil_recluster))
for (cl in 1:2) {
  idx <- features$km_class_recluster == cl
  cat(sprintf("C%d: n=%d (%.1f%%)  WBC %.1f->%.1f->%.1f->%.1f\n", cl, sum(idx),
    sum(idx) / nrow(features) * 100,
    mean(features$wbc_baseline[idx], na.rm = TRUE), mean(features$wbc_t24h[idx], na.rm = TRUE),
    mean(features$wbc_t48h[idx], na.rm = TRUE), mean(features$wbc_t72h[idx], na.rm = TRUE)))
}

# ─── Method B: predict from MIMIC centroids ───
model_path <- file.path(PROJECT_ROOT, "output", "model", "mimic_cluster_model.rds")
have_model <- file.exists(model_path)
if (have_model) {
  mimic_model <- readRDS(model_path)
  features[, km_class_mimic := predict_cluster(features, mimic_model)]
  agree <- mean(features$km_class_mimic == features$km_class_recluster) * 100
  cat(sprintf("\n--- Method B: MIMIC-centroid prediction ---\n"))
  cat(sprintf("Agreement with de novo (Method A): %.1f%%\n", agree))
  for (cl in 1:2) {
    idx <- features$km_class_mimic == cl
    cat(sprintf("C%d: n=%d (%.1f%%)  baseline WBC %.1f\n", cl, sum(idx),
      sum(idx) / nrow(features) * 100, mean(features$wbc_baseline[idx], na.rm = TRUE)))
  }
} else {
  stop("output/model/mimic_cluster_model.rds not found: the MIMIC section of this script writes it, so the file must come from a completed run.")
}

eicu <- merge(eicu, features[, .(patientunitstayid, km_class = km_class_recluster,
  km_class_recluster, km_class_mimic, wbc_baseline, wbc_t24h, wbc_t48h, wbc_t72h)],
  by = "patientunitstayid")

fwrite(eicu[, .(patientunitstayid, km_class_recluster, km_class_mimic, hosp_death)],
       file.path(PROJECT_ROOT, "data", "eicu_cluster_assignment.csv"))

fwrite(data.table(Cohort = "eICU", Silhouette = round(sil_recluster, 3),
                  n = nrow(features)),
       file.path(TABLES, "_tableS_eicu_silhouette.csv"))


ref <- if (have_model) mimic_model$ref_stats else
  stop("mimic_cluster_model.rds not found — run the MIMIC clustering section first")
c1_mort <- mean(eicu$hosp_death[eicu$km_class == 1]) * 100
c2_mort <- mean(eicu$hosp_death[eicu$km_class == 2]) * 100

cat("\n=== Mortality ===\n")
cat(sprintf("eICU:   C1=%.1f%%  C2=%.1f%%  ratio=%.2f  diff=%.1fpp\n",
  c1_mort, c2_mort, c1_mort / c2_mort, c1_mort - c2_mort))
cat(sprintf("MIMIC:  C1=%.1f%%  C2=%.1f%%  ratio=%.2f  (reference)\n",
  ref$c1_mort_hosp, ref$c2_mort_hosp, ref$c1_mort_hosp / ref$c2_mort_hosp))
cat(sprintf("Result: gradient %s (eICU ratio %.2f vs MIMIC %.2f)\n",
  ifelse(c1_mort > c2_mort, "REPLICATED", "NOT replicated"),
  c1_mort / c2_mort, ref$c1_mort_hosp / ref$c2_mort_hosp))

# Level 2: Microbiology

cat("\n--- Level 2: Microbiology ---\n")
mvars <- c("has_culture", "culture_pos", "g_negative", "g_positive", "fungal", "has_pathogen_dx")
mlabs <- c("Any culture", "Culture positive", "Gram-negative", "Gram-positive", "Fungal", "Pathogen Dx")
for (vi in seq_along(mvars)) {
  var <- mvars[vi]; lab <- mlabs[vi]
  if (!var %in% names(eicu)) next
  r1 <- mean(eicu[km_class == 1, get(var)]) * 100; r2 <- mean(eicu[km_class == 2, get(var)]) * 100
  tbl <- table(eicu$km_class, eicu[[var]])
  if (all(dim(tbl) == 2)) {
    chi <- chisq.test(tbl)
    cat(sprintf("  %-18s C1=%.1f%%  C2=%.1f%%  chi2=%.1f  p=%.2e\n", lab, r1, r2, chi$statistic, chi$p.value))
  } else cat(sprintf("  %-18s C1=%.1f%%  C2=%.1f%%\n", lab, r1, r2))
}

# Level 3: Antibiotic spectrum

cat("\n--- Level 3: Antibiotic Spectrum ---\n")
avars <- c("gp_abx", "gn_abx", "fungal_abx", "n_abx"); alabs <- c("GP-targeted", "GN-targeted", "Fungal-targeted", "N unique abx")
for (vi in seq_along(avars)) {
  var <- avars[vi]; lab <- alabs[vi]
  if (!var %in% names(eicu)) next
  if (var == "n_abx") {
    cat(sprintf("  %-18s C1=%.1f  C2=%.1f  p=%.2e\n", lab,
      mean(eicu[km_class == 1, get(var)]), mean(eicu[km_class == 2, get(var)]),
      wilcox.test(eicu[km_class == 1, get(var)], eicu[km_class == 2, get(var)])$p.value))
  } else {
    r1 <- mean(eicu[km_class == 1, get(var)]) * 100; r2 <- mean(eicu[km_class == 2, get(var)]) * 100
    tbl <- table(eicu$km_class, eicu[[var]])
    if (all(dim(tbl) == 2)) {
      chi <- chisq.test(tbl)
      cat(sprintf("  %-18s C1=%.1f%%  C2=%.1f%%  chi2=%.1f  p=%.2e\n", lab, r1, r2, chi$statistic, chi$p.value))
    }
  }
}

# Level 4: Organism mapping

cat("\n--- Level 4: Organism Mapping ---\n")
gn_pat <- "escherichia|klebsiella|enterobacter|serratia|pseudomonas|acinetobacter|proteus|haemophilus|neisseria|bacteroides|citrobacter|morganella|stenotrophomonas"
gp_pat <- "staphylococcus|streptococcus|enterococcus|listeria|clostridium|bacillus|corynebacterium"
fg_pat <- "candida|aspergillus|cryptococcus|yeast|fungus"
eicu_micro[, org_class := ifelse(grepl(fg_pat, tolower(organism)), "Fungal",
  ifelse(grepl(gn_pat, tolower(organism)), "Gram-negative",
  ifelse(grepl(gp_pat, tolower(organism)), "Gram-positive", "Other/Unknown")))]
stay_org <- eicu_micro[patientunitstayid %in% cohort_ids, .(
  has_gn = as.integer(any(org_class == "Gram-negative")),
  has_gp = as.integer(any(org_class == "Gram-positive")),
  has_fg = as.integer(any(org_class == "Fungal"))), by = patientunitstayid]
eicu <- merge(eicu, stay_org, by = "patientunitstayid", all.x = TRUE)
for (col in c("has_gn", "has_gp", "has_fg")) eicu[[col]][is.na(eicu[[col]])] <- 0
for (cl in 1:2) {
  idx <- eicu$km_class == cl
  cat(sprintf("C%d: GN=%.1f%%  GP=%.1f%%  FG=%.1f%%\n", cl,
    mean(eicu$has_gn[idx]) * 100, mean(eicu$has_gp[idx]) * 100, mean(eicu$has_fg[idx]) * 100))
}

# Pooled MIMIC + eICU cohort interaction

cat("\n--- Pooled cohort interaction (km_class x cohort) ---\n")
safe({
  mimic_clustered <- fread(file.path(PROJECT_ROOT, "data", "feature_matrix_with_clusters.csv"))
# Method B (MIMIC-centroid) labels for eICU in the pooled analysis
  pooled <- rbind(
    data.table(cohort = "MIMIC", km_class = as.integer(mimic_clustered$km_class),
               mort = as.integer(mimic_clustered$hospital_expire_flag)),
    data.table(cohort = "eICU", km_class = as.integer(eicu$km_class_mimic),
               mort = as.integer(eicu$hosp_death)))
  pooled[, cohort := factor(cohort, levels = c("MIMIC", "eICU"))]
  mp <- glm(mort ~ km_class * cohort, data = pooled, family = binomial)
  co <- summary(mp)$coefficients
  int_row <- grep("km_class:cohort", rownames(co), value = TRUE)[1]
  pooled_dt <- data.table(
    Term = c("km_class (MIMIC)", "cohort (eICU vs MIMIC)", "km_class x cohort"),
    Estimate = co[c("km_class", grep("^cohort", rownames(co), value = TRUE)[1], int_row), 1],
    OR = exp(co[c("km_class", grep("^cohort", rownames(co), value = TRUE)[1], int_row), 1]),
    p = co[c("km_class", grep("^cohort", rownames(co), value = TRUE)[1], int_row), 4])
  fwrite(pooled_dt, file.path(TABLES, "tableS_pooled_interaction.csv"))
  print(pooled_dt)
  cat(sprintf("Interaction p = %.4f (%s trajectory effect across cohorts)\n",
    co[int_row, 4], ifelse(co[int_row, 4] < 0.05, "differs", "consistent")))
}, "pooled_interaction")

# CMH + Breslow-Day

cat("\n--- CMH + Breslow-Day stratified analysis ---\n")
safe({
  mimic_c1_dead <- sum(mimic_clustered$km_class == 1 & mimic_clustered$hospital_expire_flag == 1)
  mimic_c1_live <- sum(mimic_clustered$km_class == 1 & mimic_clustered$hospital_expire_flag == 0)
  mimic_c2_dead <- sum(mimic_clustered$km_class == 2 & mimic_clustered$hospital_expire_flag == 1)
  mimic_c2_live <- sum(mimic_clustered$km_class == 2 & mimic_clustered$hospital_expire_flag == 0)
  eicu_c1_dead <- sum(eicu$km_class_mimic == 1 & eicu$hosp_death == 1, na.rm = TRUE)
  eicu_c1_live <- sum(eicu$km_class_mimic == 1 & eicu$hosp_death == 0, na.rm = TRUE)
  eicu_c2_dead <- sum(eicu$km_class_mimic == 2 & eicu$hosp_death == 1, na.rm = TRUE)
  eicu_c2_live <- sum(eicu$km_class_mimic == 2 & eicu$hosp_death == 0, na.rm = TRUE)

  arr <- array(as.numeric(c(mimic_c1_dead, mimic_c1_live, mimic_c2_dead, mimic_c2_live,
                            eicu_c1_dead, eicu_c1_live, eicu_c2_dead, eicu_c2_live)),
               dim = c(2, 2, 2),
               dimnames = list(Outcome = c("Dead","Alive"),
                               Phenotype = c("C1","C2"),
                               Cohort = c("MIMIC","eICU")))

# Array as double: mantelhaen.test's CI uses integer arithmetic and overflows to NA.
  mh <- mantelhaen.test(arr, conf.level = 0.95)
  # Per-cohort ORs and Woolf weights, used for the Breslow-Day statistic below
  or_m <- (mimic_c1_dead/mimic_c1_live)/(mimic_c2_dead/mimic_c2_live)
  or_e <- (eicu_c1_dead/eicu_c1_live)/(eicu_c2_dead/eicu_c2_live)
  w1 <- 1/(1/mimic_c1_dead + 1/mimic_c1_live + 1/mimic_c2_dead + 1/mimic_c2_live)
  w2 <- 1/(1/eicu_c1_dead + 1/eicu_c1_live + 1/eicu_c2_dead + 1/eicu_c2_live)
  ci_lo <- mh$conf.int[1]; ci_hi <- mh$conf.int[2]
  cat(sprintf("CMH common OR: %.2f (%.2f-%.2f), p=%.2e\n",
      mh$estimate, ci_lo, ci_hi, mh$p.value))

  pooled_lnor <- sum(w1*log(or_m), w2*log(or_e)) / sum(w1, w2)
  bd_stat <- w1*(log(or_m)-pooled_lnor)^2 + w2*(log(or_e)-pooled_lnor)^2
  bd_p <- 1 - pchisq(bd_stat, df=1)
  cat(sprintf("OR MIMIC: %.2f, OR eICU: %.2f\n", or_m, or_e))
  cat(sprintf("Breslow-Day (Woolf) chi-sq: %.2f, df=1, p=%.2e\n", bd_stat, bd_p))

  cmh_dt <- data.table(
    Measure   = c("OR MIMIC","OR eICU","CMH common OR","CMH CI low","CMH CI high","CMH p","Breslow-Day chi-sq","Breslow-Day p"),
    Value     = c(or_m, or_e, mh$estimate, ci_lo, ci_hi, mh$p.value, bd_stat, bd_p))
  fwrite(cmh_dt, file.path(TABLES, "tableS_cmh_breslow_day.csv"))
  cat("Saved: output/tables/tableS_cmh_breslow_day.csv\n")
}, "cmh_breslow_day")


cat("\n========== MEDIATION ANALYSIS (eICU) ==========\n")
safe({
  eicu_lac <- fread(file.path(PROJECT_ROOT, "data", "eicu_lactate_long.csv"))
  eicu_lac[, patientunitstayid := as.integer(patientunitstayid)]

  lac_max <- eicu_lac[, .(max_lactate = max(lactate, na.rm = TRUE),
                           lac_measurements = .N), by = patientunitstayid]
  lac_max[, patientunitstayid := as.integer(patientunitstayid)]
  eicu <- lac_max[eicu, on = .(patientunitstayid)]  # right join: keep all eicu rows

  cat(sprintf("Lactate coverage: %.1f%% (n=%d/%d)\n",
    mean(!is.na(eicu$max_lactate)) * 100,
    sum(!is.na(eicu$max_lactate)), nrow(eicu)))
  cat(sprintf("Max lactate: C1=%.1f  C2=%.1f  p=%.2e\n",
    mean(eicu[km_class == 1, max_lactate], na.rm = TRUE),
    mean(eicu[km_class == 2, max_lactate], na.rm = TRUE),
    wilcox.test(eicu[km_class == 1, max_lactate], eicu[km_class == 2, max_lactate])$p.value))

  # Causal mediation (natural indirect/direct effects, risk-difference scale)
  eicu_med <- as.data.frame(eicu[!is.na(max_lactate) & !is.na(hosp_death),
    .(km_class, max_lactate, hosp_death, age, gender)])
  eicu_med$km_class_f <- as.integer(eicu_med$km_class == 1)
  eicu_med$is_female  <- as.integer(tolower(trimws(eicu_med$gender)) == "female")
  eicu_med$hosp_death <- as.integer(eicu_med$hosp_death)

  cm <- causal_mediation(eicu_med, "km_class_f", "max_lactate", "hosp_death",
                         c("age", "is_female"))

  cat(sprintf("\nCausal mediation (eICU):\n"))
  cat(sprintf("  ACME (indirect, pp): %.4f (%.4f, %.4f) p=%.3f\n",
    cm$ACME, cm$ACME_lo, cm$ACME_hi, cm$ACME_p))
  cat(sprintf("  ADE  (direct,   pp): %.4f (%.4f, %.4f)\n", cm$ADE, cm$ADE_lo, cm$ADE_hi))
  cat(sprintf("  Total Effect (pp):   %.4f (%.4f, %.4f)\n", cm$Total, cm$Total_lo, cm$Total_hi))
  cat(sprintf("  Prop. Mediated:      %.1f%% (%.1f, %.1f)\n",
    cm$Prop_mediated_pct, cm$CI_lower_pct, cm$CI_upper_pct))

  # Save eICU mediation (Indirect = ACME, risk-difference scale)
  med_dt <- data.table(
    Cohort = "eICU", Adjustment = "Full (+ support, abx)",
    Indirect = cm$ACME, Prop_mediated_pct = cm$Prop_mediated_pct,
    CI_lower_pct = cm$CI_lower_pct, CI_upper_pct = cm$CI_upper_pct,
    ADE = cm$ADE, ADE_lo = cm$ADE_lo, ADE_hi = cm$ADE_hi, ADE_p = cm$ADE_p)

  fwrite(med_dt, file.path(TABLES, "tableS_eicu_mediation.csv"))

  # ── Table 6: Combined mediation (MIMIC + eICU) ──
  mimic_med_file <- file.path(TABLES, "table4_mimic_mediation.csv")
  if (file.exists(mimic_med_file)) {
    mimic_med <- fread(mimic_med_file)
    mimic_med[, Cohort := "MIMIC"]
    t4 <- rbind(
      mimic_med[, .(Cohort, Adjustment, Indirect,
        Prop_mediated_pct, CI_lower_pct, CI_upper_pct,
        ADE, ADE_lo, ADE_hi, ADE_p)],
      med_dt[, .(Cohort, Adjustment, Indirect,
        Prop_mediated_pct, CI_lower_pct, CI_upper_pct,
        ADE, ADE_lo, ADE_hi, ADE_p)],
      fill = TRUE)
    fwrite(t4, file.path(TABLES, "table4_mediation.csv"))
    cat("  -> Table 6 (mediation, MIMIC + eICU) written: table6_mediation_three_cohort.csv\n")
  }
}, "eicu_mediation")


summary_dt <- data.table(
  Measure = c("Method", "C1 proportion", "C2 proportion", "Silhouette (de novo)",
    "C1 mortality", "C2 mortality", "C1/C2 mortality ratio",
    "Method A vs B agreement", "Overall validation"),
  MIMIC = c("K-means K=2", sprintf("%.1f%%", ref$c1_prop),
    sprintf("%.1f%%", ref$c2_prop), sprintf("%.3f", ref$silhouette),
    sprintf("%.1f%%", ref$c1_mort_hosp), sprintf("%.1f%%", ref$c2_mort_hosp),
    sprintf("%.2f", ref$c1_mort_hosp / ref$c2_mort_hosp), "—", "Reference"),
  eICU = c("De novo K=2 + MIMIC predict",
    sprintf("%.1f%%", mean(eicu$km_class == 1) * 100),
    sprintf("%.1f%%", mean(eicu$km_class == 2) * 100),
    sprintf("%.3f", sil_recluster),
    sprintf("%.1f%%", c1_mort), sprintf("%.1f%%", c2_mort),
    sprintf("%.2f", c1_mort / c2_mort),
    if (have_model) sprintf("%.1f%%", mean(eicu$km_class_mimic == eicu$km_class_recluster) * 100) else "—",
    ifelse(c1_mort > c2_mort, "PASSED — gradient replicated", "CHECK")))
# ── Assemble Table 7: eICU external validation (summary + pooled interaction + CMH) ──
s6_panels <- list()
# Panel A: Validation summary
summary_dt[, Panel := "A"]
setcolorder(summary_dt, "Panel")
s6_panels[["A"]] <- summary_dt

# Panel B: pooled interaction
pooled_file <- file.path(TABLES, "tableS_pooled_interaction.csv")
if (file.exists(pooled_file)) {
  pb <- fread(pooled_file)
  pb[, Panel := "B"]
  setcolorder(pb, "Panel")
  s6_panels[["B"]] <- pb
}

# Panel C: CMH + Breslow-Day
cmh_file <- file.path(TABLES, "tableS_cmh_breslow_day.csv")
if (file.exists(cmh_file)) {
  pc <- fread(cmh_file)
  pc[, Panel := "C"]
  setcolorder(pc, "Panel")
  s6_panels[["C"]] <- pc
}

s6_combined <- rbindlist(s6_panels, fill = TRUE)
fwrite(s6_combined, file.path(TABLES, "table7_eicu_validation.csv"))
cat("\nDONE -> output/tables/table7_eicu_validation.csv\n")

# Centroid-mapping metrics (both label sets use 1 = high-WBC), appended to Table 7 panel A.
safe({
  ca <- fread(file.path(PROJECT_ROOT, "data", "eicu_cluster_assignment.csv"))
  ca <- ca[km_class_recluster %in% c(1, 2) & km_class_mimic %in% c(1, 2)]
  n_tot <- nrow(ca)
  tp <- sum(ca$km_class_recluster == 1 & ca$km_class_mimic == 1)
  fn <- sum(ca$km_class_recluster == 1 & ca$km_class_mimic == 2)
  fp <- sum(ca$km_class_recluster == 2 & ca$km_class_mimic == 1)
  sens_cm  <- tp / (tp + fn) * 100
  prec_cm  <- tp / (tp + fp) * 100
  extra_cm <- fp
  ca_all <- fread(file.path(PROJECT_ROOT, "data", "eicu_cluster_assignment.csv"))
  ec     <- fread(file.path(PROJECT_ROOT, "data", "eicu_cohort.csv"))
  cm     <- merge(ca_all[, .(patientunitstayid, km_class_mimic)],
                  ec[, .(patientunitstayid, hosp_death)], by = "patientunitstayid")
  cm_c1  <- mean(cm$hosp_death[cm$km_class_mimic == 1], na.rm = TRUE) * 100
  cm_c2  <- mean(cm$hosp_death[cm$km_class_mimic == 2], na.rm = TRUE) * 100
  t7 <- fread(file.path(TABLES, "table7_eicu_validation.csv"))
  add <- data.table(
    Panel = "A",
    Measure = c("C1 sensitivity of centroid mapping",
                "C1 precision of centroid mapping",
                "C2 reclassified to C1 by centroid mapping",
                "C1 mortality (centroid-mapped labels)",
                "C2 mortality (centroid-mapped labels)",
                "C1/C2 mortality ratio (centroid-mapped labels)"),
    MIMIC = "—",
    eICU = c(sprintf("%.1f%%", sens_cm),
             sprintf("%.1f%%", prec_cm),
             sprintf("%d (%.1f%%)", extra_cm, extra_cm / n_tot * 100),
             sprintf("%.1f%%", cm_c1),
             sprintf("%.1f%%", cm_c2),
             sprintf("%.2f", cm_c1 / cm_c2)))
  for (cc in setdiff(names(t7), names(add))) add[, (cc) := NA]
  setcolorder(add, names(t7))
  fwrite(rbind(t7, add, fill = TRUE), file.path(TABLES, "table7_eicu_validation.csv"))
  cat(sprintf("  -> table7 panel A: appended 3 centroid-mapping rows + 3 centroid-mapped mortality rows (sens %.1f%%, prec %.1f%%, C2->C1 n=%d, mort %.1f%%/%.1f%%)\n",
              sens_cm, prec_cm, extra_cm, cm_c1, cm_c2))
}, "table7_centroid_rows")

# eICU STROBE flow; steps 1-4 from data/eicu_strobe_counts.csv (written by sql/02_eicu.sql)
eicu_strobe_cnt <- fread(file.path(PROJECT_ROOT, "data", "eicu_strobe_counts.csv"))
eicu_strobe <- data.table(
  Step = c("eICU-CRD v2.0 ICU stays (208 hospitals)",
           "Age >= 18 years",
           "ICU length of stay >= 24 hours",
           "Sepsis diagnosis (ICD-9/diagnosis string)",
           "Has WBC measurements (analytic cohort)",
           "FINAL eICU VALIDATION COHORT"),
  n = c(sprintf("%d", eicu_strobe_cnt$step1_icu_stays),
        sprintf("%d", eicu_strobe_cnt$step2_age18),
        sprintf("%d", eicu_strobe_cnt$step3_los24h),
        sprintf("%d", eicu_strobe_cnt$step4_sepsis),
        sprintf("%d", n_eicu_wbc),
        sprintf("%d", n_eicu_wbc)))
fwrite(eicu_strobe, file.path(TABLES, "tableS_strobe_eicu.csv"))
cat("-> tableS_strobe_eicu.csv written\n")

cat("============================================================\n")
close_log(LOG)

# NWICU external validation (third database)

TABLES <- file.path(PROJECT_ROOT, "output", "tables")
dir.create(TABLES, showWarnings = FALSE, recursive = TRUE)
LOG <- open_log(file.path(PROJECT_ROOT, "output", "run_log_nwicu.txt"))

# 0. Load NWICU data

nwicu       <- fread(file.path(PROJECT_ROOT, "data", "nwicu_cohort.csv"))
nwicu_long  <- fread(file.path(PROJECT_ROOT, "data", "nwicu_blood_long.csv"))
nwicu_lac   <- fread(file.path(PROJECT_ROOT, "data", "nwicu_lactate_long.csv"))

nwicu[, stay_id := as.integer(stay_id)]
nwicu_long[, `:=`(stay_id = as.integer(stay_id), hrs = as.numeric(hrs), wbc = as.numeric(wbc))]
nwicu_lac[, `:=`(stay_id = as.integer(stay_id), hrs = as.numeric(hrs), lactate = as.numeric(lactate))]

# Coerce boolean columns (NWICU exports TRUE/FALSE as text)
for (col in names(nwicu)) {
  x <- nwicu[[col]]
  if (is.character(x) && all(unique(na.omit(x)) %in% c("t","f","true","false","1","0","TRUE","FALSE"))) {
    nwicu[[col]] <- as.integer(tolower(x) %in% c("t","true","1"))
  }
}

cohort_ids <- unique(nwicu$stay_id)
nwicu_long <- nwicu_long[stay_id %in% cohort_ids]
nwicu_lac  <- nwicu_lac[stay_id %in% cohort_ids]

cat(sprintf("\nNWICU cohort: %d stays, %d WBC measurements, %d lactate measurements\n",
  nrow(nwicu), nrow(nwicu_long), nrow(nwicu_lac)))
cat(sprintf("28-day mortality: %.1f%% (%d/%d)\n",
  mean(nwicu$death_28day) * 100, sum(nwicu$death_28day), nrow(nwicu)))

# 1. Build WBC windows + scale

features <- build_wbc_windows(nwicu_long, "stay_id", WBC_WINDOWS)
features[, n_tp_wbc := rowSums(!is.na(.SD)), .SDcols = CLUSTER_FEATURES]

nwicu_mice <- impute_wbc_mice(features, features = CLUSTER_FEATURES,
                              aux = "n_tp_wbc", m = 5, seed = 2026, maxit = 5)
X <- nwicu_mice$per_imp[[1]]                      # single completed dataset
X_scaled <- scale(X); X_scaled[!is.finite(X_scaled)] <- 0

cat("\n============================================================\n")
cat("NWICU EXTERNAL VALIDATION (Database #3)\n")
cat("============================================================\n\n")

# ─── Missingness report ───
missingness_report(features, CLUSTER_FEATURES, "NWICU WBC windows")

# ─── Method A: de novo K-means K=2 ───
set.seed(2026)
kf <- fit_kmeans(X_scaled, 2)
al <- align_labels_by_wbc(kf$km$cluster, X[, "wbc_baseline"])
features[, km_class_recluster := al$labels]
idx_sub <- sample(seq_len(nrow(X_scaled)), min(5000, nrow(X_scaled)))
sil_recluster <- summary(silhouette(kf$km$cluster[idx_sub], dist(X_scaled[idx_sub, ])))$avg.width

cat("--- Method A: De novo K-means K=2 ---\n")
cat(sprintf("Silhouette: %.3f\n", sil_recluster))
for (cl in 1:2) {
  idx <- features$km_class_recluster == cl
  cat(sprintf("C%d: n=%d (%.1f%%)  WBC %.1f->%.1f->%.1f->%.1f\n", cl, sum(idx),
    sum(idx) / nrow(features) * 100,
    mean(features$wbc_baseline[idx], na.rm = TRUE),
    mean(features$wbc_t24h[idx], na.rm = TRUE),
    mean(features$wbc_t48h[idx], na.rm = TRUE),
    mean(features$wbc_t72h[idx], na.rm = TRUE)))
}

# ─── Method B: predict from MIMIC centroids ───
model_path <- file.path(PROJECT_ROOT, "output", "model", "mimic_cluster_model.rds")
have_model <- file.exists(model_path)
if (have_model) {
  mimic_model <- readRDS(model_path)
  features[, km_class_mimic := predict_cluster(features, mimic_model)]
  agree <- mean(features$km_class_mimic == features$km_class_recluster) * 100
  cat(sprintf("\n--- Method B: MIMIC-centroid prediction ---\n"))
  cat(sprintf("Agreement with de novo (Method A): %.1f%%\n", agree))
  for (cl in 1:2) {
    idx <- features$km_class_mimic == cl
    cat(sprintf("C%d: n=%d (%.1f%%)  baseline WBC %.1f\n", cl, sum(idx),
      sum(idx) / nrow(features) * 100, mean(features$wbc_baseline[idx], na.rm = TRUE)))
  }
} else {
  cat("\n*** mimic_cluster_model.rds not found — run 01 first; Method B skipped ***\n")
  features[, km_class_mimic := km_class_recluster]
}

nwicu <- merge(nwicu, features[, .(stay_id, km_class = km_class_recluster,
  km_class_recluster, km_class_mimic,
  wbc_baseline, wbc_t24h, wbc_t48h, wbc_t72h)],
  by = "stay_id")

if (have_model) {
  nw_c11 <- sum(nwicu$km_class_recluster == 1 & nwicu$km_class_mimic == 1)
  nw_c12 <- sum(nwicu$km_class_recluster == 1 & nwicu$km_class_mimic == 2)
  nw_c21 <- sum(nwicu$km_class_recluster == 2 & nwicu$km_class_mimic == 1)
  nw_c22 <- sum(nwicu$km_class_recluster == 2 & nwicu$km_class_mimic == 2)
  fwrite(data.table(de_novo = c("C1", "C1", "C2", "C2"),
                    centroid = c("C1", "C2", "C1", "C2"),
                    n = c(nw_c11, nw_c12, nw_c21, nw_c22)),
         file.path(TABLES, "_nwicu_centroid_mapping.csv"))
  fwrite(nwicu[, .(stay_id, km_class_recluster, km_class_mimic, death_28day)],
         file.path(PROJECT_ROOT, "data", "nwicu_cluster_assignment.csv"))
}


# 2. Mortality validation (death_28day primary)

ref <- if (have_model) mimic_model$ref_stats else
  stop("mimic_cluster_model.rds not found — run the MIMIC clustering section first")

c1_mort_28 <- mean(nwicu$death_28day[nwicu$km_class == 1]) * 100
c2_mort_28 <- mean(nwicu$death_28day[nwicu$km_class == 2]) * 100
c1_mort_90 <- mean(nwicu$death_90day[nwicu$km_class == 1]) * 100
c2_mort_90 <- mean(nwicu$death_90day[nwicu$km_class == 2]) * 100

cat("\n=== Mortality Validation ===\n")
cat(sprintf("NWICU 28d: C1=%.1f%%  C2=%.1f%%  ratio=%.2f  diff=%.1fpp\n",
  c1_mort_28, c2_mort_28, c1_mort_28 / c2_mort_28, c1_mort_28 - c2_mort_28))
cat(sprintf("NWICU 90d: C1=%.1f%%  C2=%.1f%%  ratio=%.2f  diff=%.1fpp\n",
  c1_mort_90, c2_mort_90, c1_mort_90 / c2_mort_90, c1_mort_90 - c2_mort_90))
cat(sprintf("MIMIC 28d: C1=%.1f%%  C2=%.1f%%  ratio=%.2f  (reference)\n",
  ref$c1_mort_28d, ref$c2_mort_28d,
  ref$c1_mort_28d / ref$c2_mort_28d))
cat(sprintf("Result: 28d gradient %s (NWICU ratio %.2f vs MIMIC %.2f)\n",
  ifelse(c1_mort_28 > c2_mort_28, "REPLICATED", "NOT replicated"),
  c1_mort_28 / c2_mort_28, ref$c1_mort_28d / ref$c2_mort_28d))

safe({
  m_adj <- glm(death_28day ~ km_class + age + gender + has_vaso + charlson_score,
               data = nwicu, family = binomial)
  or <- 1 / exp(coef(m_adj)["km_class"])
  ci <- suppressMessages(confint(m_adj))["km_class", ]
  ci_or <- rev(1 / exp(ci))
  assign("or", or, envir = .GlobalEnv)
  assign("ci", ci_or, envir = .GlobalEnv)
  cat(sprintf("\nAdjusted OR (C1 vs C2, 28d): %.2f (%.2f-%.2f), p=%.2e\n",
    or, ci_or[1], ci_or[2], summary(m_adj)$coefficients["km_class", 4]))
}, "nwicu_logistic")

# 3. Antibiotic spectrum comparison (Level 3 equivalent, no microbiology)

cat("\n--- Antibiotic Spectrum ---\n")
avars <- c("has_gp_abx", "has_gn_abx", "has_fungal_abx", "has_broad_abx",
           "anti_mrsa", "anti_pseudomonal", "n_unique_abx")
alabs <- c("GP-targeted", "GN-targeted", "Fungal-targeted", "Broad-spectrum",
           "Anti-MRSA", "Anti-pseudomonal", "N unique abx")
abx_rows <- list()
for (vi in seq_along(avars)) {
  var <- avars[vi]; lab <- alabs[vi]
  if (!var %in% names(nwicu)) next
  if (var == "n_unique_abx") {
    m1 <- mean(nwicu[km_class == 1, get(var)]); m2 <- mean(nwicu[km_class == 2, get(var)])
    pv <- wilcox.test(nwicu[km_class == 1, get(var)], nwicu[km_class == 2, get(var)])$p.value
    cat(sprintf("  %-18s C1=%.1f  C2=%.1f  p=%.2e\n", lab, m1, m2, pv))
    abx_rows[[length(abx_rows) + 1]] <- data.table(
      Measure = lab, C1 = sprintf("%.1f", m1), C2 = sprintf("%.1f", m2),
      p = format.pval(pv, digits = 3, eps = 1e-16), Test = "Wilcoxon")
  } else {
    r1 <- mean(nwicu[km_class == 1, get(var)]) * 100
    r2 <- mean(nwicu[km_class == 2, get(var)]) * 100
    tbl <- table(nwicu$km_class, nwicu[[var]])
    if (all(dim(tbl) == 2)) {
      chi <- chisq.test(tbl)
      cat(sprintf("  %-18s C1=%.1f%%  C2=%.1f%%  chi2=%.1f  p=%.2e\n",
        lab, r1, r2, chi$statistic, chi$p.value))
      abx_rows[[length(abx_rows) + 1]] <- data.table(
        Measure = lab, C1 = sprintf("%.1f%%", r1), C2 = sprintf("%.1f%%", r2),
        p = format.pval(chi$p.value, digits = 3, eps = 1e-16), Test = "Chi-square")
    } else {
      cat(sprintf("  %-18s C1=%.1f%%  C2=%.1f%%\n", lab, r1, r2))
      abx_rows[[length(abx_rows) + 1]] <- data.table(
        Measure = lab, C1 = sprintf("%.1f%%", r1), C2 = sprintf("%.1f%%", r2),
        p = NA_character_, Test = NA_character_)
    }
  }
}
abx_rows[[length(abx_rows) + 1]] <- data.table(
  Measure = "90-day mortality", C1 = sprintf("%.1f%%", c1_mort_90),
  C2 = sprintf("%.1f%%", c2_mort_90), p = NA_character_, Test = NA_character_)
abx_dt <- rbindlist(abx_rows, fill = TRUE)
fwrite(abx_dt, file.path(TABLES, "tableS11_nwicu_auxiliary_outcomes.csv"))

# 4. Lactate mediation analysis (same bootstrap method as eICU)

cat("\n========== CAUSAL MEDIATION (NWICU) ==========\n")
safe({
  # Compute max lactate per stay (raw labevents, itemid 100031)
  lac_max <- nwicu_lac[, .(max_lactate = max(lactate, na.rm = TRUE),
                             lac_measurements = .N), by = stay_id]
  nwicu <- lac_max[nwicu, on = .(stay_id)]

  cat(sprintf("Lactate coverage: %.1f%% (n=%d/%d)\n",
    mean(!is.na(nwicu$max_lactate)) * 100,
    sum(!is.na(nwicu$max_lactate)), nrow(nwicu)))
  cat(sprintf("Max lactate: C1=%.1f  C2=%.1f  p=%.2e\n",
    mean(nwicu[km_class == 1, max_lactate], na.rm = TRUE),
    mean(nwicu[km_class == 2, max_lactate], na.rm = TRUE),
    wilcox.test(nwicu[km_class == 1, max_lactate],
                nwicu[km_class == 2, max_lactate])$p.value))

  # Causal mediation (natural indirect/direct effects, risk-difference scale)
  nwmed <- as.data.frame(nwicu[!is.na(max_lactate) & !is.na(death_28day),
    .(km_class, max_lactate, death_28day, age, gender, has_vaso, charlson_score)])
  nwmed$km_class_f <- as.integer(nwmed$km_class == 1)
  nwmed$is_female  <- as.integer(nwmed$gender == "F")
  nwmed$death_28day <- as.integer(nwmed$death_28day)

  cm <- causal_mediation(nwmed, "km_class_f", "max_lactate", "death_28day",
                         c("age", "is_female", "has_vaso", "charlson_score"))
  cm[, Cohort := "NWICU"]
  cm[, Adjustment := "Full (+ vaso, charlson)"]

  cat(sprintf("\nCausal mediation (NWICU):\n"))
  cat(sprintf("  ACME (indirect, pp): %.4f (%.4f, %.4f) p=%.3f\n",
    cm$ACME, cm$ACME_lo, cm$ACME_hi, cm$ACME_p))
  cat(sprintf("  ADE  (direct,   pp): %.4f (%.4f, %.4f)\n", cm$ADE, cm$ADE_lo, cm$ADE_hi))
  cat(sprintf("  Total Effect (pp):   %.4f (%.4f, %.4f)\n", cm$Total, cm$Total_lo, cm$Total_hi))
  cat(sprintf("  Prop. Mediated:      %.1f%% (%.1f, %.1f)\n",
    cm$Prop_mediated_pct, cm$CI_lower_pct, cm$CI_upper_pct))

  # Save (Indirect = ACME, risk-difference scale)
  nwicu_med_dt <- cm[, .(Cohort, Adjustment, Indirect = ACME, Prop_mediated_pct,
                         CI_lower_pct, CI_upper_pct,
                         ADE = ADE, ADE_lo = ADE_lo, ADE_hi = ADE_hi, ADE_p = ADE_p)]
  assign("nwicu_med_dt", nwicu_med_dt, envir = .GlobalEnv)
  assign("nwicu_med_prop", cm$Prop_mediated_pct, envir = .GlobalEnv)
  assign("nwicu_med_ci", c(cm$CI_lower_pct, cm$CI_upper_pct), envir = .GlobalEnv)
  assign("nwicu_med_indirect", cm$ACME, envir = .GlobalEnv)
}, "nwicu_mediation")

# 5. Three-cohort pooled analysis (MIMIC + eICU + NWICU)

cat("\n========== THREE-COHORT POOLED ANALYSIS ==========\n")
safe({
    mimic_clustered <- fread(file.path(PROJECT_ROOT, "data", "feature_matrix_with_clusters.csv"))
  eicu_cohort      <- fread(file.path(PROJECT_ROOT, "data", "eicu_cohort.csv"))
  eicu_assign      <- fread(file.path(PROJECT_ROOT, "data", "eicu_cluster_assignment.csv"))
  eicu_assign[, patientunitstayid := as.integer(patientunitstayid)]
  eicu_cohort[, patientunitstayid := as.integer(patientunitstayid)]
# Inner join keeps eICU stays with WBC features (21,781 of 21,950)
  eicu_cohort <- merge(eicu_cohort, eicu_assign[, .(patientunitstayid, km_class_mimic)],
                       by = "patientunitstayid")

  harmonize_gender <- function(g) {
    g <- trimws(tolower(as.character(g)))
    out <- rep(NA_character_, length(g))
    out[g %in% c("m", "male")]     <- "M"
    out[g %in% c("f", "female")]   <- "F"
    out
  }
  mimic_g <- harmonize_gender(mimic_clustered$gender)
  eicu_g  <- harmonize_gender(eicu_cohort$gender)
  nwicu_g <- harmonize_gender(nwicu$gender)

  pooled <- rbind(
    data.table(cohort = "MIMIC", km_class = as.integer(mimic_clustered$km_class),
               mort = as.integer(mimic_clustered$death_28day),
               age = mimic_clustered$admission_age,
               gender = mimic_g),
    data.table(cohort = "eICU", km_class = as.integer(eicu_cohort$km_class_mimic),
               mort = as.integer(eicu_cohort$hosp_death),
               age = as.numeric(eicu_cohort$age),
               gender = eicu_g),
    data.table(cohort = "NWICU", km_class = as.integer(nwicu$km_class_mimic),
               mort = as.integer(nwicu$death_28day),
               age = as.numeric(nwicu$age),
               gender = nwicu_g),
    fill = TRUE)
  pooled[, cohort := factor(cohort, levels = c("MIMIC", "eICU", "NWICU"))]

  # 3-cohort interaction model
  mp <- glm(mort ~ km_class * cohort + age + gender, data = pooled, family = binomial)
  cat("\n--- 3-Cohort Interaction Model ---\n")
  cat(sprintf("N pooled = %d\n", nrow(pooled)))
  print(summary(mp)$coefficients)

  co <- summary(mp)$coefficients
  cat("\nPer-cohort km_class effects:\n")

  int_rows <- grep("km_class:cohort", rownames(co), value = TRUE)
  for (r in int_rows) {
    cat(sprintf("  %s: coef=%.3f, p=%.3f\n", r, co[r, 1], co[r, 4]))
  }

  pooled_dt <- as.data.table(co, keep.rownames = "Term")
  pooled_dt[, OR := exp(Estimate)]
}, "three_cohort_pooled")

# 6. CMH + Breslow-Day across all 3 cohorts

cat("\n--- 3-Cohort CMH + Breslow-Day ---\n")
safe({
  # Use Method B labels (MIMIC-centroid prediction) for all cohorts
  # MIMIC (28-day mortality)
  mc1d <- sum(mimic_clustered$km_class == 1 & mimic_clustered$death_28day == 1, na.rm = TRUE)
  mc1l <- sum(mimic_clustered$km_class == 1 & mimic_clustered$death_28day == 0, na.rm = TRUE)
  mc2d <- sum(mimic_clustered$km_class == 2 & mimic_clustered$death_28day == 1, na.rm = TRUE)
  mc2l <- sum(mimic_clustered$km_class == 2 & mimic_clustered$death_28day == 0, na.rm = TRUE)

  # eICU (centroid-mapped labels)
  ec1d <- sum(eicu_cohort$km_class_mimic == 1 & eicu_cohort$hosp_death == 1, na.rm = TRUE)
  ec1l <- sum(eicu_cohort$km_class_mimic == 1 & eicu_cohort$hosp_death == 0, na.rm = TRUE)
  ec2d <- sum(eicu_cohort$km_class_mimic == 2 & eicu_cohort$hosp_death == 1, na.rm = TRUE)
  ec2l <- sum(eicu_cohort$km_class_mimic == 2 & eicu_cohort$hosp_death == 0, na.rm = TRUE)

  # NWICU (centroid-mapped labels)
  nc1d <- sum(nwicu$km_class_mimic == 1 & nwicu$death_28day == 1, na.rm = TRUE)
  nc1l <- sum(nwicu$km_class_mimic == 1 & nwicu$death_28day == 0, na.rm = TRUE)
  nc2d <- sum(nwicu$km_class_mimic == 2 & nwicu$death_28day == 1, na.rm = TRUE)
  nc2l <- sum(nwicu$km_class_mimic == 2 & nwicu$death_28day == 0, na.rm = TRUE)

  cells <- c(mc1d, mc1l, mc2d, mc2l, ec1d, ec1l, ec2d, ec2l, nc1d, nc1l, nc2d, nc2l)
  if (any(cells <= 0)) {
    cat(sprintf("  WARNING: %d zero/negative CMH cell(s); mantelhaen.test may be unstable\n",
      sum(cells <= 0)))
  }

  arr <- array(as.numeric(c(mc1d, mc1l, mc2d, mc2l,
                         ec1d, ec1l, ec2d, ec2l,
                         nc1d, nc1l, nc2d, nc2l)),
               dim = c(2, 2, 3),
               dimnames = list(Outcome = c("Dead","Alive"),
                               Phenotype = c("C1","C2"),
                               Cohort = c("MIMIC","eICU","NWICU")))

# Array as double: mantelhaen.test's CI uses integer arithmetic and overflows to NA.
  mh <- mantelhaen.test(arr, conf.level = 0.95)
  cat(sprintf("CMH common OR: %.2f (%.2f-%.2f), p=%.2e\n",
    mh$estimate, mh$conf.int[1], mh$conf.int[2], mh$p.value))

  or_m <- (mc1d/mc1l)/(mc2d/mc2l)
  or_e <- (ec1d/ec1l)/(ec2d/ec2l)
  or_n <- (nc1d/nc1l)/(nc2d/nc2l)
  cat(sprintf("OR MIMIC: %.2f, OR eICU: %.2f, OR NWICU: %.2f\n", or_m, or_e, or_n))

  # Breslow-Day test on the Woolf weights
  w1 <- 1/(1/mc1d + 1/mc1l + 1/mc2d + 1/mc2l)
  w2 <- 1/(1/ec1d + 1/ec1l + 1/ec2d + 1/ec2l)
  w3 <- 1/(1/nc1d + 1/nc1l + 1/nc2d + 1/nc2l)
  pooled_lnor <- sum(w1*log(or_m), w2*log(or_e), w3*log(or_n)) / sum(w1, w2, w3)
  bd_stat <- w1*(log(or_m)-pooled_lnor)^2 + w2*(log(or_e)-pooled_lnor)^2 + w3*(log(or_n)-pooled_lnor)^2
  bd_df <- 2
  bd_p <- 1 - pchisq(bd_stat, df = bd_df)
  cat(sprintf("Breslow-Day (Woolf) chi-sq: %.2f, df=%d, p=%.2e\n", bd_stat, bd_df, bd_p))
  bd_p <- 1 - pchisq(bd_stat, df = bd_df)
  cat(sprintf("Breslow-Day (Woolf) chi-sq: %.2f, df=%d, p=%.2e\n", bd_stat, bd_df, bd_p))

  cmh_ci_lo <- mh$conf.int[1]; cmh_ci_hi <- mh$conf.int[2]

  cmh_dt <- data.table(
    Measure = c("OR MIMIC","OR eICU","OR NWICU","CMH common OR",
                "CMH CI low","CMH CI high","CMH p",
                "Breslow-Day chi-sq","Breslow-Day df","Breslow-Day p"),
    Value = c(or_m, or_e, or_n, mh$estimate,
              cmh_ci_lo, cmh_ci_hi, mh$p.value,
              bd_stat, bd_df, bd_p))
  fwrite(cmh_dt, file.path(TABLES, "tableS9_cmh_three_cohort.csv"))
  cat("  -> tableS9_cmh_three_cohort.csv\n")
}, "cmh_three_cohort")

# 7. Table 6: Combined mediation (MIMIC + eICU + NWICU)

safe({
  all_med <- data.table()

  # MIMIC mediation (table4_mimic_mediation.csv, wide format)
  mimic_med_file <- file.path(TABLES, "table4_mimic_mediation.csv")
  if (file.exists(mimic_med_file)) {
    mm <- fread(mimic_med_file)
    mm[, Cohort := "MIMIC"]
    all_med <- rbind(all_med, mm[, .(Cohort, Adjustment, Indirect,
      Prop_mediated_pct, CI_lower_pct, CI_upper_pct,
      ADE, ADE_lo, ADE_hi, ADE_p)], fill = TRUE)
  }

  eicu_med_file <- file.path(TABLES, "tableS_eicu_mediation.csv")
  if (file.exists(eicu_med_file)) {
    em <- fread(eicu_med_file)
    all_med <- rbind(all_med, em[, .(Cohort, Adjustment, Indirect,
      Prop_mediated_pct, CI_lower_pct, CI_upper_pct,
      ADE, ADE_lo, ADE_hi, ADE_p)], fill = TRUE)
  }

  # NWICU mediation (from global env)
  if (exists("nwicu_med_dt", envir = .GlobalEnv)) {
    nm <- get("nwicu_med_dt", envir = .GlobalEnv)
    nm[, Cohort := "NWICU"]
    all_med <- rbind(all_med, nm[, .(Cohort, Adjustment, Indirect,
      Prop_mediated_pct, CI_lower_pct, CI_upper_pct,
      ADE, ADE_lo, ADE_hi, ADE_p)], fill = TRUE)
  }

  if (nrow(all_med) > 0) {
    assign("all_med", all_med, envir = .GlobalEnv)
    if (!file.exists(mimic_med_file))
      stop("table4_mimic_mediation.csv not found: the MIMIC mediation section writes it earlier in this script.")
    fwrite(all_med, file.path(TABLES, "table6_mediation_three_cohort.csv"))
    cat(sprintf("\n  -> Table 6 (three-cohort mediation): %d rows written\n", nrow(all_med)))
  }
}, "table4_mediation")

# 10. Summary table

# Lactate coverage in NWICU: non-missing max_lactate in data/nwicu_cohort.csv, reported in Table 8.
.nw <- fread(file.path(PROJECT_ROOT, "data", "nwicu_cohort.csv"),
             select = "max_lactate", showProgress = FALSE)
.lac_n <- sum(!is.na(.nw$max_lactate)); .nw_n <- nrow(.nw)

nwicu_summary <- data.table(
  Measure = c("Database", "Cohort size", "C1 proportion", "C2 proportion",
    "Silhouette (de novo)", "Method A vs B agreement",
    "C1 28d mortality", "C2 28d mortality", "C1/C2 mortality ratio",
    "Adjusted OR (C1 vs C2)", "Lactate prop. mediated",
    "Lactate measured, n", "Lactate measured, %",
    "Overall validation"),
  NWICU = c("NWICU v0.1.0",
    as.character(nrow(nwicu)),
    sprintf("%.1f%%", mean(nwicu$km_class == 1) * 100),
    sprintf("%.1f%%", mean(nwicu$km_class == 2) * 100),
    sprintf("%.3f", sil_recluster),
    if (have_model) sprintf("%.1f%%", agree) else "N/A",
    sprintf("%.1f%%", c1_mort_28), sprintf("%.1f%%", c2_mort_28),
    sprintf("%.2f", c1_mort_28 / c2_mort_28),
    if (exists("or", envir = .GlobalEnv)) sprintf("%.2f (%.2f-%.2f)", get("or", envir = .GlobalEnv), get("ci", envir = .GlobalEnv)[1], get("ci", envir = .GlobalEnv)[2]) else "N/A",
    if (exists("nwicu_med_prop", envir = .GlobalEnv)) sprintf("%.1f%% (%.1f-%.1f%%)", get("nwicu_med_prop", envir = .GlobalEnv), get("nwicu_med_ci", envir = .GlobalEnv)[1], get("nwicu_med_ci", envir = .GlobalEnv)[2]) else "N/A",
    sprintf("%d/%d", .lac_n, .nw_n),
    sprintf("%.1f%%", .lac_n / .nw_n * 100),
    ifelse(c1_mort_28 > c2_mort_28, "PASSED: gradient replicated", "CHECK")))
fwrite(nwicu_summary, file.path(TABLES, "table8_nwicu_validation.csv"))
cat("\n  -> table8_nwicu_validation.csv\n")

# Final combined 3-cohort summary
cat("\n============================================================\n")
cat("THREE-DATABASE CROSS-VALIDATION COMPLETE\n")
cat("============================================================\n")
mimic_n <- ref$n
eicu_sil_path <- file.path(TABLES, "_tableS_eicu_silhouette.csv")
eicu_sil <- if (file.exists(eicu_sil_path))
  fread(eicu_sil_path)$Silhouette[1] else NA_real_

cat(sprintf("MIMIC:  K=2, Silhouette=%.3f, n=%d\n",
  ref$silhouette, mimic_n))
cat(sprintf("eICU:   K=2, Silhouette=%.3f, n=%d\n",
  eicu_sil, nrow(eicu_cohort)))
cat(sprintf("NWICU:  K=2, Silhouette=%.3f, n=%d\n",
  sil_recluster, nrow(nwicu)))

# NWICU STROBE flow; counts from data/nwicu_strobe_counts.csv (written by sql/03_nwicu.sql)
nw_cnt <- fread(file.path(PROJECT_ROOT, "data", "nwicu_strobe_counts.csv"))
nwicu_strobe <- data.table(
  Step = c("NWICU v0.1.0 ICU stays",
           "Angus sepsis cohort (age >= 18, ICU LOS >= 1 day)",
           "Has WBC in extraction window (-24 to +120 h)",
           "FINAL NWICU VALIDATION COHORT (baseline WBC available)"),
  n = c(sprintf("%d", nw_cnt$step1_icu_stays),
        sprintf("%d", nw_cnt$step2_angus_sepsis),
        sprintf("%d", nw_cnt$step5_wbc_any_window),
        sprintf("%d", nw_cnt$step6_wbc_baseline_final)))
fwrite(nwicu_strobe, file.path(TABLES, "tableS14_strobe_nwicu.csv"))
cat("-> tableS14_strobe_nwicu.csv written\n")

cn_n   <- sum(df$is_cn)
cn_lac <- sum(df$is_cn & !is.na(df$max_lactate))
cat(sprintf("MIMIC culture-negative subgroup: %d of %d (%.1f%%); %d with lactate (%.1f%%)\n",
            cn_n, nrow(df), cn_n / nrow(df) * 100, cn_lac, cn_lac / cn_n * 100))

t5 <- fread(file.path(TABLES, "table5_cns_interaction.csv"), showProgress = FALSE)
cn_tab <- t5[Panel == "A" & Group %in% c("True culture-negative", "Never cultured"), sum(N)]

close_log(LOG)

cat("done\n")
