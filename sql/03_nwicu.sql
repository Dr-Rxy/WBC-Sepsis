-- NWICU (external validation cohort 2)


-- PHASE 1: Anchor cohort + exclusion flags
DROP TABLE IF EXISTS tmp_nw_cohort;
CREATE TEMP TABLE tmp_nw_cohort AS
SELECT
    sc.stay_id, sc.hadm_id, sc.subject_id,
    sc.admission_age, sc.gender,
    sc.first_careunit, sc.intime, sc.outtime, sc.icu_los_days,
    sc.deathtime, sc.dischtime, sc.intime AS admittime,
    sc.admission_year,
    oc.death_28day, oc.death_90day,
    COALESCE(vaso.has_vaso, FALSE) AS has_vaso,
    COALESCE(lab.max_creatinine, 0) AS max_creatinine,
    COALESCE(ch.charlson_score, 0) AS charlson_score
FROM derived.sepsis_cohort sc
LEFT JOIN derived.outcomes     oc  ON sc.stay_id = oc.stay_id
LEFT JOIN derived.vasopressor  vaso ON sc.stay_id = vaso.stay_id
LEFT JOIN derived.sofa_labs    lab ON sc.stay_id = lab.stay_id
LEFT JOIN derived.charlson     ch  ON sc.stay_id = ch.stay_id
WHERE sc.admission_age >= 18
  AND sc.icu_los_days >= 1;
CREATE INDEX ON tmp_nw_cohort (stay_id);
CREATE INDEX ON tmp_nw_cohort (hadm_id);

-- PHASE 2: Time zero = first IV antibiotic per stay
DROP TABLE IF EXISTS tmp_nw_timezero;
CREATE TEMP TABLE tmp_nw_timezero AS
SELECT
    c.stay_id, c.hadm_id,
    MIN(abx.starttime) AS time_zero,
    COUNT(DISTINCT abx.drug) AS n_abx_drugs
FROM tmp_nw_cohort c
JOIN derived.antibiotic abx ON c.stay_id = abx.stay_id
WHERE abx.starttime IS NOT NULL
GROUP BY c.stay_id, c.hadm_id;
CREATE INDEX ON tmp_nw_timezero (stay_id);

-- PHASE 3: WBC longitudinal — align to time_zero
DROP TABLE IF EXISTS tmp_nw_wbc_long;
CREATE TEMP TABLE tmp_nw_wbc_long AS
SELECT
    tz.stay_id,
    tz.hadm_id,
    tz.time_zero,
    le.charttime::timestamp AS charttime,
    le.valuenum::numeric AS wbc,
    EXTRACT(EPOCH FROM (le.charttime::timestamp - tz.time_zero)) / 3600.0 AS hrs_from_tz
FROM tmp_nw_timezero tz
JOIN nw_hosp.labevents le ON tz.hadm_id = le.hadm_id::integer
WHERE le.itemid = '100016'
  AND le.valuenum IS NOT NULL
  AND le.valuenum !~ '[^0-9.]'              -- exclude non-numeric junk
  AND le.valuenum::numeric > 0              -- physiologically plausible
  AND le.valuenum::numeric < 500            -- extreme outlier cap
  AND le.charttime::timestamp
      BETWEEN tz.time_zero - INTERVAL '24 hours'
      AND tz.time_zero + INTERVAL '120 hours';  -- generous extraction window
CREATE INDEX ON tmp_nw_wbc_long (stay_id);
CREATE INDEX ON tmp_nw_wbc_long (hrs_from_tz);

-- PHASE 4: WBC window features (same 4 windows as MIMIC + eICU)
-- WBC windows (hrs from time_zero): [-8,24) [12,36) [36,60) [60,84)
DROP TABLE IF EXISTS tmp_nw_wbc_features;
CREATE TEMP TABLE tmp_nw_wbc_features AS
SELECT
    stay_id,
    COUNT(*) AS n_wbc_measurements,
    MIN(hrs_from_tz) AS first_wbc_hrs,
    MAX(hrs_from_tz) AS last_wbc_hrs,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY wbc)
        FILTER (WHERE hrs_from_tz >= -8  AND hrs_from_tz < 24)  AS wbc_baseline,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY wbc)
        FILTER (WHERE hrs_from_tz >= 12 AND hrs_from_tz < 36)  AS wbc_t24h,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY wbc)
        FILTER (WHERE hrs_from_tz >= 36 AND hrs_from_tz < 60)  AS wbc_t48h,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY wbc)
        FILTER (WHERE hrs_from_tz >= 60 AND hrs_from_tz < 84)  AS wbc_t72h
FROM tmp_nw_wbc_long
GROUP BY stay_id;
CREATE INDEX ON tmp_nw_wbc_features (stay_id);

-- PHASE 5: Lactate longitudinal (for mediation analysis)
DROP TABLE IF EXISTS tmp_nw_lactate_long;
CREATE TEMP TABLE tmp_nw_lactate_long AS
SELECT
    tz.stay_id,
    le.charttime::timestamp AS charttime,
    le.valuenum::numeric AS lactate,
    EXTRACT(EPOCH FROM (le.charttime::timestamp - tz.time_zero)) / 3600.0 AS hrs_from_tz
FROM tmp_nw_timezero tz
JOIN nw_hosp.labevents le ON tz.hadm_id = le.hadm_id::integer
WHERE le.itemid = '100031'                   -- Lactate
  AND le.valuenum IS NOT NULL
  AND le.valuenum !~ '[^0-9.]'
  AND le.valuenum::numeric > 0
  AND le.valuenum::numeric < 100             -- extreme cap
  AND le.charttime::timestamp
      BETWEEN tz.time_zero - INTERVAL '24 hours'
      AND tz.time_zero + INTERVAL '120 hours';
CREATE INDEX ON tmp_nw_lactate_long (stay_id);

-- Max lactate per stay from raw labevents (itemid 100031, mmol/L), aligned to time_zero.
DROP TABLE IF EXISTS tmp_nw_lactate_max;
CREATE TEMP TABLE tmp_nw_lactate_max AS
SELECT stay_id, MAX(lactate) AS max_lactate
FROM tmp_nw_lactate_long
GROUP BY stay_id;
CREATE INDEX ON tmp_nw_lactate_max (stay_id);

-- PHASE 6: Antibiotic spectrum features (replacing microbiology)
DROP TABLE IF EXISTS tmp_nw_abx_features;
CREATE TEMP TABLE tmp_nw_abx_features AS
SELECT
    tz.stay_id,
    COUNT(DISTINCT abx.drug) AS n_unique_abx,
    BOOL_OR(LOWER(abx.drug) ~ 'vancomycin|linezolid|daptomycin') AS has_gp_abx,
    BOOL_OR(LOWER(abx.drug) ~ 'piperacillin|cefepime|ceftazidime|meropenem|imipenem|aztreonam|ciprofloxacin|levofloxacin|gentamicin|tobramycin|amikacin') AS has_gn_abx,
    BOOL_OR(LOWER(abx.drug) ~ 'fluconazole|caspofungin|micafungin|amphotericin|voriconazole|posaconazole') AS has_fungal_abx,
    BOOL_OR(LOWER(abx.drug) ~ 'meropenem|imipenem|doripenem|cefepime|ceftazidime|piperacillin|vancomycin|linezolid|daptomycin|ciprofloxacin|levofloxacin|amikacin|tobramycin|gentamicin') AS has_broad_abx,
    BOOL_OR(LOWER(abx.drug) ~ 'vancomycin|linezolid|daptomycin') AS anti_mrsa,
    BOOL_OR(LOWER(abx.drug) ~ 'piperacillin|cefepime|ceftazidime|meropenem|imipenem|aztreonam') AS anti_pseudomonal
FROM tmp_nw_timezero tz
JOIN derived.antibiotic abx ON tz.stay_id = abx.stay_id
WHERE abx.starttime <= tz.time_zero + INTERVAL '72 hours'
GROUP BY tz.stay_id;
CREATE INDEX ON tmp_nw_abx_features (stay_id);

-- PHASE 7: Assemble final cohort CSV
\copy (SELECT c.stay_id, c.hadm_id, c.subject_id, c.admission_age AS age, c.gender, c.first_careunit, c.icu_los_days, tz.time_zero, tz.n_abx_drugs, COALESCE(wf.wbc_baseline, wf.wbc_t24h) AS wbc_baseline, wf.wbc_t24h, wf.wbc_t48h, wf.wbc_t72h, wf.n_wbc_measurements AS n_wbc_tp, wf.first_wbc_hrs, wf.last_wbc_hrs, lac.max_lactate, c.max_creatinine, c.has_vaso, c.charlson_score, COALESCE(abx.n_unique_abx, 0) AS n_unique_abx, COALESCE(abx.has_gp_abx, FALSE) AS has_gp_abx, COALESCE(abx.has_gn_abx, FALSE) AS has_gn_abx, COALESCE(abx.has_fungal_abx, FALSE) AS has_fungal_abx, COALESCE(abx.has_broad_abx, FALSE) AS has_broad_abx, COALESCE(abx.anti_mrsa, FALSE) AS anti_mrsa, COALESCE(abx.anti_pseudomonal, FALSE) AS anti_pseudomonal, c.death_28day, c.death_90day, c.deathtime, c.dischtime, c.admission_year FROM tmp_nw_cohort c JOIN tmp_nw_timezero tz ON c.stay_id = tz.stay_id LEFT JOIN tmp_nw_wbc_features wf ON c.stay_id = wf.stay_id LEFT JOIN tmp_nw_abx_features abx ON c.stay_id = abx.stay_id LEFT JOIN tmp_nw_lactate_max lac ON c.stay_id = lac.stay_id WHERE wf.wbc_baseline IS NOT NULL AND tz.time_zero IS NOT NULL ORDER BY c.stay_id) TO 'data/nwicu_cohort.csv' WITH CSV HEADER;

-- PHASE 8: Export blood_long CSV (longitudinal format for trajectory plots)
\copy (SELECT wbc.stay_id, ROUND(wbc.hrs_from_tz::numeric, 2) AS hrs, wbc.wbc FROM tmp_nw_wbc_long wbc WHERE wbc.stay_id IN ( SELECT stay_id FROM tmp_nw_wbc_features WHERE wbc_baseline IS NOT NULL ) ORDER BY wbc.stay_id, wbc.hrs_from_tz) TO 'data/nwicu_blood_long.csv' WITH CSV HEADER;

-- PHASE 9: Export lactate_long CSV (for mediation bootstrap)
\copy (SELECT lac.stay_id, ROUND(lac.hrs_from_tz::numeric, 2) AS hrs, lac.lactate FROM tmp_nw_lactate_long lac WHERE lac.stay_id IN ( SELECT stay_id FROM tmp_nw_wbc_features WHERE wbc_baseline IS NOT NULL ) ORDER BY lac.stay_id, lac.hrs_from_tz) TO 'data/nwicu_lactate_long.csv' WITH CSV HEADER;

-- Diagnostics
SELECT '─── NWICU WBC Extraction Diagnostics ───' AS stage;
SELECT 'Cohort (adult, LOS>=1)' AS step, COUNT(*) AS n FROM tmp_nw_cohort
UNION ALL SELECT 'With time zero', COUNT(*) FROM tmp_nw_timezero
UNION ALL SELECT 'With WBC features', COUNT(*) FROM tmp_nw_wbc_features
UNION ALL SELECT 'WBC baseline coverage', COUNT(*) FROM tmp_nw_wbc_features WHERE wbc_baseline IS NOT NULL
UNION ALL SELECT 'WBC t24h coverage', COUNT(*) FROM tmp_nw_wbc_features WHERE wbc_t24h IS NOT NULL
UNION ALL SELECT 'WBC t48h coverage', COUNT(*) FROM tmp_nw_wbc_features WHERE wbc_t48h IS NOT NULL
UNION ALL SELECT 'WBC t72h coverage', COUNT(*) FROM tmp_nw_wbc_features WHERE wbc_t72h IS NOT NULL;

SELECT 'WBC measurements per stay' AS metric,
    ROUND(AVG(n_wbc_measurements), 1) AS mean,
    ROUND(MIN(n_wbc_measurements)::numeric, 0) AS min,
    ROUND(PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY n_wbc_measurements)::numeric, 0) AS median,
    ROUND(MAX(n_wbc_measurements)::numeric, 0) AS max
FROM tmp_nw_wbc_features;

-- ============================================================================
-- NWICU STROBE flow: progressive cohort counts
-- ============================================================================

-- ── Step 1: all ICU stays ───────────────────────────────────────────────────
DROP TABLE IF EXISTS st_s1;
CREATE TEMP TABLE st_s1 AS SELECT COUNT(*)::int AS n FROM nw_icu.icustays;

-- ── Steps 2-3: Angus sepsis cohort, then age / LOS filter ───────────────────
DROP TABLE IF EXISTS st_s2;
CREATE TEMP TABLE st_s2 AS SELECT COUNT(*)::int AS n FROM derived.sepsis_cohort;

DROP TABLE IF EXISTS st_cohort;
CREATE TEMP TABLE st_cohort AS
SELECT sc.stay_id, sc.hadm_id
FROM derived.sepsis_cohort sc
WHERE sc.admission_age >= 18
  AND sc.icu_los_days >= 1;
CREATE INDEX ON st_cohort (stay_id);
CREATE INDEX ON st_cohort (hadm_id);

DROP TABLE IF EXISTS st_s3;
CREATE TEMP TABLE st_s3 AS SELECT COUNT(*)::int AS n FROM st_cohort;

-- ── Step 4: time zero = first IV antibiotic per stay ───────────────────────
DROP TABLE IF EXISTS st_tz;
CREATE TEMP TABLE st_tz AS
SELECT c.stay_id, MIN(abx.starttime) AS tz
FROM st_cohort c
JOIN derived.antibiotic abx ON abx.stay_id = c.stay_id
WHERE abx.starttime IS NOT NULL
GROUP BY c.stay_id;
CREATE INDEX ON st_tz (stay_id);

DROP TABLE IF EXISTS st_s4;
CREATE TEMP TABLE st_s4 AS SELECT COUNT(*)::int AS n FROM st_tz;

-- ── Single pass over the 16.7 M-row labevents table ────────────────────────
DROP TABLE IF EXISTS st_wbc_raw;
CREATE TEMP TABLE st_wbc_raw AS
SELECT le.hadm_id, le.charttime::timestamp AS charttime, le.valuenum::numeric AS wbc
FROM nw_hosp.labevents le
WHERE le.itemid = '100016'
  AND le.valuenum IS NOT NULL
  AND le.valuenum !~ '[^0-9.]'
  AND le.valuenum::numeric > 0
  AND le.valuenum::numeric < 500;
CREATE INDEX ON st_wbc_raw (hadm_id);

-- ── Step 5: WBC present anywhere in the extraction window (-24..+120 h) ────
DROP TABLE IF EXISTS st_s5;
CREATE TEMP TABLE st_s5 AS
SELECT COUNT(DISTINCT t.stay_id)::int AS n
FROM st_tz t
JOIN st_cohort c ON c.stay_id = t.stay_id
JOIN st_wbc_raw w ON w.hadm_id = c.hadm_id::text
WHERE w.charttime BETWEEN t.tz - INTERVAL '24 hours'
                      AND t.tz + INTERVAL '120 hours';

-- ── Step 6 (FINAL): WBC present in the BASELINE window [-8, +24) ───────────
DROP TABLE IF EXISTS st_s6;
CREATE TEMP TABLE st_s6 AS
SELECT COUNT(DISTINCT t.stay_id)::int AS n
FROM st_tz t
JOIN st_cohort c ON c.stay_id = t.stay_id
JOIN st_wbc_raw w ON w.hadm_id = c.hadm_id::text
WHERE w.charttime BETWEEN t.tz - INTERVAL '8 hours'
                      AND t.tz + INTERVAL '24 hours';

-- -- Funnel (written to disk for the R side to read) ----------------------------
\copy (SELECT (SELECT n FROM st_s1) AS step1_icu_stays, (SELECT n FROM st_s2) AS step2_angus_sepsis, (SELECT n FROM st_s3) AS step3_age_los, (SELECT n FROM st_s4) AS step4_time_zero, (SELECT n FROM st_s5) AS step5_wbc_any_window, (SELECT n FROM st_s6) AS step6_wbc_baseline_final) TO 'data/nwicu_strobe_counts.csv' CSV HEADER

-- QC of the NWICU in-hospital death flag: its positive rate is implausible, so the primary outcome is derived.outcomes.death_28day.

-- A. Whole source table
SELECT 'A source admissions (all)' AS step,
       COUNT(*)                                               AS n,
       SUM((hospital_expire_flag)::int)                       AS n_flag1,
       ROUND(100.0 * AVG((hospital_expire_flag)::int), 3)     AS pct
FROM nw_hosp.admissions;

-- B. The cohort used here (n=4,992)
SELECT 'B Angus cohort (derived.outcomes)' AS step,
       COUNT(*)                                                AS n,
       SUM((hospital_expire_flag)::int)                        AS n_hosp_flag1,
       ROUND(100.0 * AVG((hospital_expire_flag)::int), 3)      AS pct_hosp_flag,
       SUM((death_28day)::int)                                 AS n_death_28day,
       ROUND(100.0 * AVG((death_28day)::int), 3)               AS pct_death_28day
FROM derived.outcomes;

-- C. Re-check via derived.sepsis_cohort (expect = B)
SELECT 'C Angus cohort via sepsis_cohort' AS step,
       COUNT(*)                                                   AS n,
       SUM((a.hospital_expire_flag)::int)                         AS n_flag1,
       ROUND(100.0 * AVG((a.hospital_expire_flag)::int), 3)       AS pct
FROM derived.sepsis_cohort sc
JOIN nw_hosp.admissions a ON a.hadm_id::text = sc.hadm_id::text
WHERE sc.admission_age >= 18
  AND sc.icu_los_days >= 1;

-- D. Derived ICU stay base table
SELECT 'D derived.icustay_baseline' AS step,
       COUNT(*)                                               AS n,
       SUM((hospital_expire_flag)::int)                       AS n_flag1,
       ROUND(100.0 * AVG((hospital_expire_flag)::int), 3)     AS pct
FROM derived.icustay_baseline;

-- E. Value distribution and missingness (expect 0/1, no NAs)
SELECT 'E null & distinct check' AS step,
       COUNT(*)                                                        AS n_total,
       COUNT(hospital_expire_flag)                                     AS n_notnull,
       SUM(CASE WHEN hospital_expire_flag IS NULL THEN 1 ELSE 0 END)   AS n_null,
       COUNT(DISTINCT hospital_expire_flag::text)                      AS n_distinct,
       SUM(CASE WHEN hospital_expire_flag::text = '1' THEN 1 ELSE 0 END) AS n_eq1
FROM derived.outcomes;
