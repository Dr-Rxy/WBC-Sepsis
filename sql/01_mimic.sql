-- MIMIC-IV v3.1 (derivation cohort)


-- PHASE 0: Cleanup + indexes ------------------------------------------------
DROP TABLE IF EXISTS tmp_sepsis_cohort, tmp_time_zero, tmp_blood_aligned, tmp_blood_features, tmp_labs_agg, tmp_icu, tmp_adm, tmp_sofa, tmp_vent, tmp_rrt, tmp_vaso, tmp_abx, tmp_micro, tmp_bcx, tmp_excl, tmp_scores, tmp_cohort_merged, tmp_feature_matrix;

CREATE INDEX IF NOT EXISTS idx_icustay_detail_stay     ON mimiciv_derived.icustay_detail(stay_id);
CREATE INDEX IF NOT EXISTS idx_antibiotic_stay         ON mimiciv_derived.antibiotic(stay_id);
CREATE INDEX IF NOT EXISTS idx_sofa_stay               ON mimiciv_derived.sofa(stay_id, hr);
CREATE INDEX IF NOT EXISTS idx_microbiologyevents_hadm ON mimiciv_hosp.microbiologyevents(hadm_id);
CREATE INDEX IF NOT EXISTS idx_diagnoses_icd_subject   ON mimiciv_hosp.diagnoses_icd(subject_id);
CREATE INDEX IF NOT EXISTS idx_prescriptions_hadm      ON mimiciv_hosp.prescriptions(hadm_id);
CREATE INDEX IF NOT EXISTS idx_inputevents_stay        ON mimiciv_icu.inputevents(stay_id);

-- PHASE 1: Sepsis cohort ----------------------------------------------------
CREATE TEMP TABLE tmp_sepsis_cohort AS
SELECT
    s3.subject_id, s3.stay_id, icu.hadm_id,
    s3.suspected_infection_time, s3.antibiotic_time,
    s3.culture_time, s3.sofa_time,
    s3.sofa_score AS sofa_sepsis3,
    s3.respiration, s3.coagulation, s3.liver,
    s3.cardiovascular, s3.cns, s3.renal,
    icu.intime AS icu_intime, icu.outtime AS icu_outtime,
    icu.los    AS icu_los_days,
    icu.first_careunit, icu.last_careunit
FROM mimiciv_derived.sepsis3 s3
JOIN mimiciv_icu.icustays    icu ON s3.stay_id = icu.stay_id
WHERE s3.sepsis3 = TRUE;

CREATE INDEX ON tmp_sepsis_cohort (stay_id);
CREATE INDEX ON tmp_sepsis_cohort (hadm_id);
CREATE INDEX ON tmp_sepsis_cohort (subject_id);
CREATE INDEX ON tmp_sepsis_cohort (subject_id, hadm_id);

-- PHASE 2: Time zero (antibiotic → inputevents → prescriptions → sepsis3) ---
CREATE TEMP TABLE tmp_time_zero AS
WITH abx_derived_first AS (
    SELECT stay_id, MIN(starttime) AS tz_derived
    FROM mimiciv_derived.antibiotic
    WHERE route IN ('IV', 'IVPB', 'Intravenous', 'IVP')
      AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.stay_id = antibiotic.stay_id)
    GROUP BY stay_id
),
abx_input_first AS (
    SELECT stay_id, MIN(starttime) AS tz_input
    FROM mimiciv_icu.inputevents
    WHERE EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.stay_id = inputevents.stay_id)
      AND EXISTS (SELECT 1 FROM mimiciv_icu.d_items di
                  WHERE di.itemid = inputevents.itemid
                    AND LOWER(di.label) ~* '(vancomycin|piperacillin|ceftriaxone|cefepime|meropenem|ampicillin|aztreonam|ciprofloxacin|levofloxacin|moxifloxacin|metronidazole|clindamycin|linezolid|ertapenem|imipenem|ceftazidime|cefazolin|nafcillin|oxacillin|penicillin|gentamicin|tobramycin|amikacin|tigecycline|doxycycline|sulfamethoxazole|colistin|polymyxin)'
                    AND LOWER(di.label) !~* '(peak|trough|random)')
    GROUP BY stay_id
),
abx_rx_first AS (
    SELECT hadm_id, MIN(starttime) AS tz_rx
    FROM mimiciv_hosp.prescriptions
    WHERE LOWER(route) IN ('iv', 'ivpb', 'intravenous', 'ivp')
      AND LOWER(drug) ~* '(vancomycin|piperacillin|ceftriaxone|cefepime|meropenem|ampicillin|aztreonam|ciprofloxacin|levofloxacin|moxifloxacin|metronidazole|clindamycin|linezolid|ertapenem|imipenem|ceftazidime|cefazolin|nafcillin|oxacillin|penicillin|gentamicin|tobramycin|amikacin|tigecycline|doxycycline|sulfamethoxazole|colistin|polymyxin)'
      AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.hadm_id = prescriptions.hadm_id)
    GROUP BY hadm_id
)
SELECT
    s.stay_id,
    COALESCE(abx.tz_derived, inp.tz_input, rx.tz_rx, s.antibiotic_time) AS time_zero,
    CASE WHEN abx.tz_derived IS NOT NULL THEN 'derived'
         WHEN inp.tz_input   IS NOT NULL THEN 'inputevents'
         WHEN rx.tz_rx       IS NOT NULL THEN 'prescriptions'
         ELSE 'sepsis3_fallback' END AS tz_source
FROM tmp_sepsis_cohort s
LEFT JOIN abx_derived_first abx ON s.stay_id = abx.stay_id
LEFT JOIN abx_input_first   inp ON s.stay_id = inp.stay_id
LEFT JOIN abx_rx_first       rx ON s.hadm_id = rx.hadm_id;

CREATE INDEX ON tmp_time_zero (stay_id);

-- PHASE 3a: Blood cell long-format (WBC, ANC, ALC) -------------------------
CREATE TEMP TABLE tmp_blood_aligned AS
WITH blood_raw AS (
    SELECT
        b.subject_id, b.hadm_id, b.charttime, b.wbc,
        COALESCE(b.neutrophils_abs,
            CASE WHEN b.neutrophils IS NOT NULL AND b.wbc IS NOT NULL
                 THEN b.wbc * b.neutrophils / 100.0 END) AS anc,
        COALESCE(b.lymphocytes_abs,
            CASE WHEN b.lymphocytes IS NOT NULL AND b.wbc IS NOT NULL
                 THEN b.wbc * b.lymphocytes / 100.0 END) AS alc
    FROM mimiciv_derived.blood_differential b
    WHERE b.wbc IS NOT NULL
      AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s
                  WHERE s.subject_id = b.subject_id AND s.hadm_id = b.hadm_id)

    UNION ALL

    SELECT
        le.subject_id, le.hadm_id, le.charttime, NULL AS wbc,
        CASE WHEN le.itemid IN (51255, 53153) THEN le.valuenum END AS anc,
        CASE WHEN le.itemid IN (51254, 53152) THEN le.valuenum END AS alc
    FROM mimiciv_hosp.labevents le
    WHERE le.itemid IN (51255, 53153, 51254, 53152)
      AND le.valuenum IS NOT NULL AND le.valuenum > 0 AND le.valuenum < 500
      AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.hadm_id = le.hadm_id)

    UNION ALL

    SELECT
        lp.subject_id, lp.hadm_id, lp.charttime, NULL AS wbc,
        CASE WHEN lp.itemid IN (51256, 51258)
             THEN w.wbc * lp.valuenum / 100.0 END AS anc,
        CASE WHEN lp.itemid IN (51244, 51245)
             THEN w.wbc * lp.valuenum / 100.0 END AS alc
    FROM mimiciv_hosp.labevents lp
    JOIN mimiciv_derived.blood_differential w
        ON w.subject_id = lp.subject_id AND w.hadm_id = lp.hadm_id
       AND w.charttime BETWEEN lp.charttime - INTERVAL '4 hours'
                            AND lp.charttime + INTERVAL '4 hours'
       AND w.wbc IS NOT NULL
    WHERE lp.itemid IN (51256, 51258, 51244, 51245)
      AND lp.valuenum IS NOT NULL AND lp.valuenum > 0 AND lp.valuenum <= 100
      AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.hadm_id = lp.hadm_id)
),
blood_aligned AS (
    SELECT
        s.stay_id, b.charttime, b.wbc, b.anc, b.alc, tz.time_zero,
        EXTRACT(EPOCH FROM (b.charttime - tz.time_zero)) / 3600 AS hrs
    FROM blood_raw b
    JOIN tmp_sepsis_cohort s ON b.hadm_id = s.hadm_id AND b.subject_id = s.subject_id
    JOIN tmp_time_zero    tz ON s.stay_id = tz.stay_id
    WHERE b.charttime >= tz.time_zero - INTERVAL '8 hours'
      AND b.charttime <= tz.time_zero + INTERVAL '96 hours'
)
SELECT stay_id, charttime, hrs, wbc, anc, alc FROM blood_aligned;

CREATE INDEX ON tmp_blood_aligned (stay_id);
CREATE INDEX ON tmp_blood_aligned (stay_id, hrs);

-- PHASE 3b: Aggregate to time windows ---------------------------------------
-- V3 disjoint grid: [-8,24) [24,48) [48,72) [72,96), named by START hour.
-- Half-open intervals throughout, so the four windows partition [-8,96) and no
-- measurement enters two features. The published grid used BETWEEN (closed, so
-- boundary measurements were counted twice) with a 32-h baseline window (12-24 h
-- shared with wbc_t24h). See docs/MODIFICATIONS.md.
-- The published wbc_t96h / anc_t96h / alc_t96h window (84-108 h) is REMOVED: the
-- extraction stops at +96 h, so that window is entirely contained in wbc_t72h and
-- can no longer serve as an independent imputation auxiliary.
CREATE TEMP TABLE tmp_blood_features AS
SELECT
    stay_id,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY wbc) FILTER (WHERE hrs >= -8 AND hrs < 24) AS wbc_baseline,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY anc) FILTER (WHERE hrs >= -8 AND hrs < 24) AS anc_baseline,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY alc) FILTER (WHERE hrs >= -8 AND hrs < 24) AS alc_baseline,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY wbc) FILTER (WHERE hrs >= 24 AND hrs < 48) AS wbc_t24h,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY anc) FILTER (WHERE hrs >= 24 AND hrs < 48) AS anc_t24h,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY alc) FILTER (WHERE hrs >= 24 AND hrs < 48) AS alc_t24h,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY wbc) FILTER (WHERE hrs >= 48 AND hrs < 72) AS wbc_t48h,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY anc) FILTER (WHERE hrs >= 48 AND hrs < 72) AS anc_t48h,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY alc) FILTER (WHERE hrs >= 48 AND hrs < 72) AS alc_t48h,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY wbc) FILTER (WHERE hrs >= 72 AND hrs < 96) AS wbc_t72h,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY anc) FILTER (WHERE hrs >= 72 AND hrs < 96) AS anc_t72h,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY alc) FILTER (WHERE hrs >= 72 AND hrs < 96) AS alc_t72h,
    -- Coverage count over the closed extraction range [-8, 96]: a measurement
    -- drawn exactly at +96 h is available and must be counted. This is NOT a
    -- window, so the half-open convention does not apply to it.
    COUNT(*) FILTER (WHERE hrs >= -8 AND hrs <= 96 AND wbc IS NOT NULL) AS n_cbc_total,
    BOOL_OR(hrs >= -8 AND hrs < 24 AND wbc IS NOT NULL) AS has_baseline_cbc
FROM tmp_blood_aligned
GROUP BY stay_id;

CREATE INDEX ON tmp_blood_features (stay_id);

-- PHASE 4: Lab aggregates ---------------------------------------------------
-- Lactate and the other laboratory maxima are taken over time_zero - 24 h to
-- + 72 h. The published extraction anchored these on ICU admission instead,
-- which put them on a different clock from the four clustering windows; eICU
-- and NWICU used yet other widths. All three cohorts now use the same
-- interval on the same anchor.
CREATE TEMP TABLE tmp_labs_agg AS
WITH labs_raw AS (
    SELECT c.stay_id, le.itemid, le.valuenum
    FROM mimiciv_hosp.labevents le
    JOIN tmp_sepsis_cohort c ON le.hadm_id = c.hadm_id
    JOIN tmp_time_zero    tz ON c.stay_id  = tz.stay_id
    WHERE le.itemid IN (50813,52442,53154,50889,51237,51675,50912,52546,50885,53089,50861,50878,51265,53189)
      AND le.valuenum IS NOT NULL AND le.valuenum > 0
      AND le.charttime >= tz.time_zero - INTERVAL '24 hours'
      AND le.charttime <= tz.time_zero + INTERVAL '72 hours'
)
SELECT stay_id,
    MAX(valuenum) FILTER (WHERE itemid IN (50813,52442,53154) AND valuenum < 30)  AS max_lactate,
    MAX(valuenum) FILTER (WHERE itemid = 50889 AND valuenum < 500)                 AS max_crp,
    MAX(valuenum) FILTER (WHERE itemid IN (51237,51675) AND valuenum < 20)         AS max_inr,
    MAX(valuenum) FILTER (WHERE itemid IN (50912,52546) AND valuenum < 50)         AS max_creatinine,
    MAX(valuenum) FILTER (WHERE itemid IN (50885,53089) AND valuenum < 100)        AS max_bilirubin,
    MAX(valuenum) FILTER (WHERE itemid = 50861 AND valuenum < 10000)               AS max_alt,
    MAX(valuenum) FILTER (WHERE itemid = 50878 AND valuenum < 10000)               AS max_ast,
    MAX(valuenum) FILTER (WHERE itemid IN (51265,53189) AND valuenum < 3000)       AS max_platelet
FROM labs_raw
GROUP BY stay_id;

CREATE INDEX ON tmp_labs_agg (stay_id);

-- PHASE 5a: Dimension tables ------------------------------------------------
DROP TABLE IF EXISTS tmp_icu, tmp_adm, tmp_sofa, tmp_org, tmp_abx, tmp_micro, tmp_bcx, tmp_excl, tmp_scores;

CREATE TEMP TABLE tmp_icu AS
SELECT stay_id, los_icu, los_hospital, admission_age,
       first_icu_stay, first_hosp_stay, gender, race,
       dod, admittime, dischtime, hospital_expire_flag
FROM mimiciv_derived.icustay_detail
WHERE EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.stay_id = icustay_detail.stay_id);
CREATE INDEX ON tmp_icu (stay_id);

CREATE TEMP TABLE tmp_adm AS
SELECT a.subject_id, a.hadm_id, a.admittime AS adm_admittime,
       a.dischtime AS adm_dischtime, a.deathtime,
       a.admission_type, a.admission_location, a.discharge_location,
       a.insurance, a.language, a.marital_status, a.race AS adm_race
FROM mimiciv_hosp.admissions a
WHERE EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.hadm_id = a.hadm_id);
CREATE INDEX ON tmp_adm (hadm_id);
CREATE INDEX ON tmp_adm (subject_id);

CREATE TEMP TABLE tmp_sofa AS
SELECT stay_id,
       MAX(respiration) AS sofa_resp, MAX(coagulation) AS sofa_coag,
       MAX(liver) AS sofa_liver, MAX(cardiovascular) AS sofa_cv,
       MAX(cns) AS sofa_cns, MAX(renal) AS sofa_renal
FROM mimiciv_derived.sofa
WHERE hr <= 24
  AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.stay_id = sofa.stay_id)
GROUP BY stay_id;
CREATE INDEX ON tmp_sofa (stay_id);

CREATE TEMP TABLE tmp_vent AS
SELECT DISTINCT stay_id, TRUE AS has_vent
FROM mimiciv_derived.ventilation
WHERE ventilation_status = 'InvasiveVent'
  AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.stay_id = ventilation.stay_id);
CREATE INDEX ON tmp_vent (stay_id);

CREATE TEMP TABLE tmp_rrt AS
SELECT DISTINCT stay_id, TRUE AS has_crrt
FROM mimiciv_derived.rrt
WHERE dialysis_present = 1
  AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.stay_id = rrt.stay_id);
CREATE INDEX ON tmp_rrt (stay_id);

CREATE TEMP TABLE tmp_vaso AS
SELECT DISTINCT stay_id, TRUE AS has_vaso
FROM mimiciv_derived.vasoactive_agent
WHERE (COALESCE(dopamine,0)>0 OR COALESCE(epinephrine,0)>0
    OR COALESCE(norepinephrine,0)>0 OR COALESCE(phenylephrine,0)>0
    OR COALESCE(vasopressin,0)>0 OR COALESCE(dobutamine,0)>0
    OR COALESCE(milrinone,0)>0)
  AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.stay_id = vasoactive_agent.stay_id);
CREATE INDEX ON tmp_vaso (stay_id);

CREATE TEMP TABLE tmp_abx AS
SELECT stay_id,
       BOOL_OR(LOWER(antibiotic) ~* 'vancomycin|teicoplanin')       AS has_glycopeptide,
       BOOL_OR(LOWER(antibiotic) ~* 'meropenem|imipenem|ertapenem') AS has_carbapenem,
       BOOL_OR(LOWER(antibiotic) ~* 'ceftriaxone|cefepime|ceftazidime') AS has_cephalosporin,
       COUNT(DISTINCT antibiotic) AS abx_n_unique
FROM mimiciv_derived.antibiotic
WHERE route IN ('IV', 'IVPB', 'Intravenous', 'IVP')
  AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.stay_id = antibiotic.stay_id)
GROUP BY stay_id;
CREATE INDEX ON tmp_abx (stay_id);

CREATE TEMP TABLE tmp_micro AS
SELECT hadm_id,
       BOOL_OR(LOWER(org_name) ~* 'escherichia coli|klebsiella|enterobacter|serratia|pseudomonas|acinetobacter|proteus|haemophilus|neisseria|bacteroides') AS has_gn,
       BOOL_OR(LOWER(org_name) ~* 'staphylococcus|streptococcus|enterococcus|listeria|clostridium') AS has_gp,
       BOOL_OR(LOWER(org_name) ~* 'candida|aspergillus|cryptococcus|yeast|fungus') AS has_fungal,
       BOOL_OR(LOWER(org_name) ~* 'candida|aspergillus|cryptococcus|yeast|fungus' AND spec_type_desc ~* 'blood culture') AS has_fungal_blood,
       BOOL_OR(LOWER(org_name) ~* 'coagulase.negative|epidermidis|hominis|corynebacterium|bacillus') AS has_contaminant,
       COUNT(DISTINCT org_name) AS n_organisms
FROM mimiciv_hosp.microbiologyevents
WHERE org_name IS NOT NULL
  AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.hadm_id = microbiologyevents.hadm_id)
GROUP BY hadm_id;
CREATE INDEX ON tmp_micro (hadm_id);

CREATE TEMP TABLE tmp_bcx AS
SELECT hadm_id,
       COUNT(DISTINCT charttime) AS n_blood_cultures
FROM mimiciv_hosp.microbiologyevents
WHERE spec_type_desc ~* 'blood culture'
  AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.hadm_id = microbiologyevents.hadm_id)
GROUP BY hadm_id;
CREATE INDEX ON tmp_bcx (hadm_id);

CREATE TEMP TABLE tmp_excl AS
SELECT DISTINCT subject_id,
       BOOL_OR(di.icd_code ~ '^C8[0-9]|^C9[0-9]|^D4[5-7]') AS excl_hemato,
       BOOL_OR(di.icd_code ~ '^(042|B2[0-4])')              AS excl_hiv,
       BOOL_OR(di.icd_code ~ '^(V42|Z94|T86)')              AS excl_transplant
FROM mimiciv_hosp.diagnoses_icd di
WHERE di.icd_code ~ '^(C8[0-9]|C9[0-9]|D4[5-7]|042|B2[0-4]|V42|Z94|T86)'
  AND EXISTS (SELECT 1 FROM tmp_sepsis_cohort s WHERE s.subject_id = di.subject_id)
GROUP BY subject_id;
CREATE INDEX ON tmp_excl (subject_id);

CREATE TEMP TABLE tmp_scores AS
SELECT c.stay_id,
       COALESCE(ap.apsiii, 0) AS apsiii,
       COALESCE(ch.charlson_comorbidity_index, 0) AS charlson
FROM tmp_sepsis_cohort c
LEFT JOIN mimiciv_derived.apsiii  ap ON c.stay_id = ap.stay_id
LEFT JOIN mimiciv_derived.charlson ch ON c.hadm_id = ch.hadm_id;
CREATE INDEX ON tmp_scores (stay_id);

-- PHASE 5b: Cohort merge + filter -------------------------------------------
DROP TABLE IF EXISTS tmp_cohort_merged;
CREATE TEMP TABLE tmp_cohort_merged AS
SELECT
    a.subject_id, a.stay_id, a.hadm_id,
    i.admission_age, i.gender, i.race,
    i.los_icu, i.los_hospital, i.first_icu_stay,
    i.dod, i.admittime, i.dischtime, i.hospital_expire_flag,
    adm.admission_type, adm.admission_location,
    a.sofa_sepsis3, a.icu_intime, a.icu_outtime, a.icu_los_days,
    a.first_careunit,
    tz.time_zero, tz.tz_source,
    sc.apsiii, sc.charlson,
    COALESCE(ex.excl_hemato, FALSE)     AS excl_hemato,
    COALESCE(ex.excl_hiv, FALSE)        AS excl_hiv,
    COALESCE(ex.excl_transplant, FALSE) AS excl_transplant,
    COALESCE(abx.has_glycopeptide, FALSE)  AS has_glycopeptide,
    COALESCE(abx.has_carbapenem, FALSE)    AS has_carbapenem,
    COALESCE(abx.has_cephalosporin, FALSE) AS has_cephalosporin,
    COALESCE(abx.abx_n_unique, 0)          AS abx_n_unique,
    COALESCE(m.has_gn, FALSE)         AS has_gram_negative,
    COALESCE(m.has_gp, FALSE)         AS has_gram_positive,
    COALESCE(m.has_fungal, FALSE)       AS has_fungal,
    COALESCE(m.has_fungal_blood, FALSE) AS has_fungal_blood,
    COALESCE(m.has_contaminant, FALSE)  AS has_contaminant,
    COALESCE(m.n_organisms, 0)         AS n_unique_organisms,
    COALESCE(bcx.n_blood_cultures, 0) AS n_blood_cultures,
    COALESCE(vent.has_vent, FALSE) AS has_vent,
    COALESCE(rrt.has_crrt, FALSE) AS has_crrt,
    COALESCE(vaso.has_vaso, FALSE) AS has_vaso,
    CASE WHEN i.dod IS NOT NULL AND i.dod <= i.admittime + INTERVAL '28 days'  THEN 1 ELSE 0 END AS death_28day,
    CASE WHEN i.dod IS NOT NULL AND i.dod <= i.admittime + INTERVAL '365 days' THEN 1 ELSE 0 END AS death_1year,
    i.dischtime - i.admittime AS hosp_los_interval
FROM tmp_sepsis_cohort a
JOIN tmp_icu           i   ON a.stay_id    = i.stay_id
JOIN tmp_time_zero     tz  ON a.stay_id    = tz.stay_id
LEFT JOIN tmp_adm      adm ON a.hadm_id    = adm.hadm_id
LEFT JOIN tmp_scores   sc  ON a.stay_id    = sc.stay_id
LEFT JOIN tmp_abx      abx ON a.stay_id    = abx.stay_id
LEFT JOIN tmp_micro    m   ON a.hadm_id    = m.hadm_id
LEFT JOIN tmp_bcx      bcx ON a.hadm_id    = bcx.hadm_id
LEFT JOIN tmp_vent    vent ON a.stay_id    = vent.stay_id
LEFT JOIN tmp_rrt      rrt ON a.stay_id    = rrt.stay_id
LEFT JOIN tmp_vaso    vaso ON a.stay_id    = vaso.stay_id
LEFT JOIN tmp_excl     ex  ON a.subject_id = ex.subject_id
WHERE i.admission_age >= 18
  AND i.first_icu_stay = TRUE
  AND i.los_icu >= 1.0
  AND tz.time_zero IS NOT NULL
  AND EXTRACT(EPOCH FROM (tz.time_zero - a.suspected_infection_time)) / 3600 BETWEEN -24 AND 72
  AND NOT COALESCE(ex.excl_hemato, FALSE)
  AND NOT COALESCE(ex.excl_hiv, FALSE)
  AND NOT COALESCE(ex.excl_transplant, FALSE);
CREATE INDEX ON tmp_cohort_merged (stay_id);
CREATE INDEX ON tmp_cohort_merged (hadm_id);
CREATE INDEX ON tmp_cohort_merged (subject_id);

-- PHASE 5c: Feature matrix --------------------------------------------------
DROP TABLE IF EXISTS tmp_feature_matrix;
CREATE TEMP TABLE tmp_feature_matrix AS
SELECT cm.*,
       bf.wbc_baseline, bf.anc_baseline, bf.alc_baseline,
       bf.wbc_t24h, bf.anc_t24h, bf.alc_t24h,
       bf.wbc_t48h, bf.anc_t48h, bf.alc_t48h,
       bf.wbc_t72h, bf.anc_t72h, bf.alc_t72h,
       bf.n_cbc_total, bf.has_baseline_cbc,
       wbc_t24h - wbc_baseline AS wbc_delta_24h,
       wbc_t48h - wbc_baseline AS wbc_delta_48h,
       wbc_t72h - wbc_baseline AS wbc_delta_72h,
       CASE WHEN wbc_baseline > 0 THEN (wbc_t24h-wbc_baseline)/wbc_baseline*100 END AS wbc_pct_24h,
       CASE WHEN wbc_baseline > 0 THEN (wbc_t72h-wbc_baseline)/wbc_baseline*100 END AS wbc_pct_72h,
       CASE WHEN wbc_baseline IS NOT NULL AND wbc_t72h IS NOT NULL THEN (wbc_t72h-wbc_baseline)/72.0 END AS wbc_slope,
       anc_t24h - anc_baseline AS anc_delta_24h,
       CASE WHEN anc_baseline IS NOT NULL AND anc_t72h IS NOT NULL THEN (anc_t72h-anc_baseline)/72.0 END AS anc_slope,
       alc_t24h - alc_baseline AS alc_delta_24h,
       CASE WHEN alc_baseline IS NOT NULL AND alc_t72h IS NOT NULL THEN (alc_t72h-alc_baseline)/72.0 END AS alc_slope,
       CASE WHEN wbc_t24h > 12 AND wbc_t48h > 12 AND wbc_t72h > 12 THEN 1 ELSE 0 END AS persistent_high_wbc,
       CASE WHEN wbc_t24h < 4  AND wbc_t48h < 4  AND wbc_t72h < 4  THEN 1 ELSE 0 END AS persistent_low_wbc,
       CASE WHEN alc_t24h < 1.0 AND alc_t48h < 1.0 AND alc_t72h < 1.0 THEN 1 ELSE 0 END AS persistent_lymphopenia,
       (CASE WHEN wbc_baseline IS NOT NULL THEN 1 ELSE 0 END +
        CASE WHEN wbc_t24h IS NOT NULL THEN 1 ELSE 0 END +
        CASE WHEN wbc_t48h IS NOT NULL THEN 1 ELSE 0 END +
        CASE WHEN wbc_t72h IS NOT NULL THEN 1 ELSE 0 END) AS n_tp_wbc
FROM tmp_cohort_merged cm
JOIN tmp_blood_features bf ON cm.stay_id = bf.stay_id
WHERE bf.n_cbc_total >= 3
  AND bf.has_baseline_cbc = TRUE;
CREATE INDEX ON tmp_feature_matrix (stay_id);

-- PHASE 5d: Final export ----------------------------------------------------
\copy (SELECT f.subject_id, f.stay_id, f.hadm_id, f.admission_age, f.gender, f.race, f.los_icu, f.los_hospital, f.first_icu_stay, f.admission_type, f.admission_location, f.icu_los_days, f.first_careunit, f.sofa_sepsis3, f.apsiii, f.charlson, f.time_zero, f.tz_source, f.has_glycopeptide, f.has_carbapenem, f.has_cephalosporin, f.abx_n_unique, f.has_gram_negative, f.has_gram_positive, f.has_fungal, f.has_fungal_blood, f.has_contaminant, f.n_unique_organisms, f.n_blood_cultures, f.has_vent, f.has_crrt, f.has_vaso, f.hospital_expire_flag, f.death_28day, f.death_1year, f.hosp_los_interval, f.wbc_baseline, f.wbc_t24h, f.wbc_t48h, f.wbc_t72h, f.anc_baseline, f.anc_t24h, f.anc_t48h, f.anc_t72h, f.alc_baseline, f.alc_t24h, f.alc_t48h, f.alc_t72h, f.wbc_delta_24h, f.wbc_delta_48h, f.wbc_delta_72h, f.wbc_pct_24h, f.wbc_pct_72h, f.wbc_slope, f.anc_delta_24h, f.anc_slope, f.alc_delta_24h, f.alc_slope, f.persistent_high_wbc, f.persistent_low_wbc, f.persistent_lymphopenia, f.n_tp_wbc, f.n_cbc_total, f.has_baseline_cbc, COALESCE(sf.sofa_resp, 0) AS sofa_respiration, COALESCE(sf.sofa_coag, 0) AS sofa_coagulation, COALESCE(sf.sofa_liver, 0) AS sofa_liver, COALESCE(sf.sofa_cv, 0) AS sofa_cardiovascular, COALESCE(sf.sofa_cns, 0) AS sofa_cns, COALESCE(sf.sofa_renal, 0) AS sofa_renal, la.max_lactate, la.max_crp, la.max_inr, la.max_creatinine, la.max_bilirubin, la.max_alt, la.max_ast, la.max_platelet FROM tmp_feature_matrix f LEFT JOIN tmp_sofa sf ON f.stay_id = sf.stay_id LEFT JOIN tmp_labs_agg la ON f.stay_id = la.stay_id ORDER BY f.subject_id, f.stay_id ) TO 'data/feature_matrix.csv' WITH CSV HEADER;

-- Export blood_long.csv -----------------------------------------------------
\copy (SELECT stay_id, hrs, wbc FROM tmp_blood_aligned WHERE wbc IS NOT NULL ORDER BY stay_id, hrs) TO 'data/blood_long.csv' WITH CSV HEADER

-- STROBE flow diagram: progressive cohort counts
\copy (WITH sepsis AS ( SELECT s3.subject_id, s3.stay_id, icu.hadm_id, s3.suspected_infection_time, s3.antibiotic_time FROM mimiciv_derived.sepsis3 s3 JOIN mimiciv_icu.icustays icu ON s3.stay_id = icu.stay_id WHERE s3.sepsis3 = TRUE ), detail AS ( SELECT s.*, d.admission_age, d.first_icu_stay, d.los_icu FROM sepsis s JOIN mimiciv_derived.icustay_detail d ON s.stay_id = d.stay_id ), tz AS ( SELECT s.stay_id, COALESCE( (SELECT MIN(a.starttime) FROM mimiciv_derived.antibiotic a WHERE a.stay_id = s.stay_id AND a.route IN ('IV','IVPB','Intravenous','IVP')), (SELECT MIN(ie.starttime) FROM mimiciv_icu.inputevents ie JOIN mimiciv_icu.d_items di ON di.itemid = ie.itemid WHERE ie.stay_id = s.stay_id AND LOWER(di.label) ~* '(vancomycin|piperacillin|ceftriaxone|cefepime|meropenem|ampicillin|aztreonam|ciprofloxacin|levofloxacin|moxifloxacin|metronidazole|clindamycin|linezolid|ertapenem|imipenem|ceftazidime|cefazolin|nafcillin|oxacillin|penicillin|gentamicin|tobramycin|amikacin|tigecycline|doxycycline|sulfamethoxazole|colistin|polymyxin)' AND LOWER(di.label) !~* '(peak|trough|random)'), (SELECT MIN(p.starttime) FROM mimiciv_hosp.prescriptions p WHERE p.hadm_id = s.hadm_id AND LOWER(p.route) IN ('iv','ivpb','intravenous','ivp') AND LOWER(p.drug) ~* '(vancomycin|piperacillin|ceftriaxone|cefepime|meropenem|ampicillin|aztreonam|ciprofloxacin|levofloxacin|moxifloxacin|metronidazole|clindamycin|linezolid|ertapenem|imipenem|ceftazidime|cefazolin|nafcillin|oxacillin|penicillin|gentamicin|tobramycin|amikacin|tigecycline|doxycycline|sulfamethoxazole|colistin|polymyxin)'), s.antibiotic_time ) AS time_zero FROM sepsis s ), excl AS ( SELECT di.subject_id, BOOL_OR(di.icd_code ~ '^C8[0-9]|^C9[0-9]|^D4[5-7]') AS hemato, BOOL_OR(di.icd_code ~ '^(042|B2[0-4])') AS hiv, BOOL_OR(di.icd_code ~ '^(V42|Z94|T86)') AS transplant FROM mimiciv_hosp.diagnoses_icd di WHERE di.icd_code ~ '^(C8[0-9]|C9[0-9]|D4[5-7]|042|B2[0-4]|V42|Z94|T86)' GROUP BY di.subject_id ), stage5 AS ( SELECT d.subject_id, d.hadm_id, d.stay_id, tz.time_zero, COALESCE(e.hemato,FALSE) AS hemato, COALESCE(e.hiv,FALSE) AS hiv, COALESCE(e.transplant,FALSE) AS transplant FROM detail d JOIN tz ON d.stay_id = tz.stay_id LEFT JOIN excl e ON d.subject_id = e.subject_id WHERE d.admission_age >= 18 AND d.first_icu_stay = TRUE AND d.los_icu >= 1.0 AND tz.time_zero IS NOT NULL AND EXTRACT(EPOCH FROM (tz.time_zero - d.suspected_infection_time)) / 3600 BETWEEN -24 AND 72 AND NOT COALESCE(e.hemato,FALSE) AND NOT COALESCE(e.hiv,FALSE) AND NOT COALESCE(e.transplant,FALSE) ) SELECT (SELECT COUNT(*) FROM mimiciv_icu.icustays) AS step1_icu_admissions, (SELECT COUNT(*) FROM sepsis) AS step2_sepsis3, (SELECT COUNT(*) FROM detail WHERE admission_age >= 18 AND first_icu_stay = TRUE AND los_icu >= 1.0) AS step3_age_first_los, (SELECT COUNT(*) FROM detail d JOIN tz ON d.stay_id = tz.stay_id WHERE d.admission_age >= 18 AND d.first_icu_stay = TRUE AND d.los_icu >= 1.0 AND tz.time_zero IS NOT NULL AND EXTRACT(EPOCH FROM (tz.time_zero - d.suspected_infection_time)) / 3600 BETWEEN -24 AND 72) AS step4_time_zero, (SELECT COUNT(*) FROM stage5) AS step5_excluded, (SELECT COUNT(*) FROM tmp_feature_matrix) AS step6_wbc_final) TO 'data/mimic_strobe_counts.csv' CSV HEADER
