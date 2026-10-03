-- eICU-CRD (external validation cohort 1)


-- PHASE 0: Indexes + sepsis cohort ------------------------------------------
CREATE INDEX IF NOT EXISTS idx_medication_stay ON eicu_icu.medication(patientunitstayid);
CREATE INDEX IF NOT EXISTS idx_microlab_stay   ON eicu_icu.microlab(patientunitstayid);

DROP TABLE IF EXISTS eicu_sepsis_cohort;
CREATE TEMP TABLE eicu_sepsis_cohort AS
SELECT DISTINCT
    p.patientunitstayid,
    p.unitdischargeoffset,
    p.hospitaldischargestatus,
    p.unittype,
    CASE WHEN p.age = '> 89' THEN 90
         WHEN p.age = '' OR p.age IS NULL THEN NULL
         ELSE CAST(p.age AS INTEGER) END AS age,
    p.gender,
    CASE WHEN p.hospitaldischargestatus = 'Expired' THEN 1 ELSE 0 END AS hosp_death
FROM eicu_icu.patient p
WHERE EXISTS (
    SELECT 1 FROM eicu_icu.diagnosis d
    WHERE d.patientunitstayid = p.patientunitstayid
      AND (LOWER(d.diagnosisstring) LIKE '%sepsis%'
           OR LOWER(d.diagnosisstring) LIKE '%septic%'
           OR LOWER(d.diagnosisstring) LIKE '%bacteremia%'
           OR d.icd9code LIKE '038%'
           OR d.icd9code LIKE '995.91%'
           OR d.icd9code LIKE '995.92%')
)
AND CASE WHEN p.age = '> 89' THEN 90
         WHEN p.age = '' OR p.age IS NULL THEN 0
         ELSE CAST(p.age AS INTEGER) END >= 18
AND p.unitdischargeoffset >= 1440;  -- ICU stay >= 24h

CREATE INDEX ON eicu_sepsis_cohort (patientunitstayid);

-- PHASE 1: Time zero --------------------------------------------------------
DROP TABLE IF EXISTS eicu_time_zero;
-- Time zero = first qualifying IV antibiotic offset; stays without one are anchored at ICU admission (offset 0).
CREATE TEMP TABLE eicu_time_zero AS
SELECT
    c.patientunitstayid,
    COALESCE(
        MIN(m.drugstartoffset) FILTER (
            WHERE LOWER(m.drugname) ~* '(vancomycin|piperacillin|ceftriaxone|cefepime|meropenem|ampicillin|aztreonam|ciprofloxacin|levofloxacin|moxifloxacin|metronidazole|clindamycin|linezolid|ertapenem|imipenem|ceftazidime|cefazolin|gentamicin|tobramycin|amikacin)'
              AND LOWER(m.routeadmin) IN ('iv','iv piggyback','intravenous','iv push')
              AND m.drugstartoffset >= -480
              AND m.drugstartoffset <= 1440
        ),
        0  -- fallback: ICU admission (offset 0)
    ) AS time_zero_min
FROM eicu_sepsis_cohort c
LEFT JOIN eicu_icu.medication m ON c.patientunitstayid = m.patientunitstayid
GROUP BY c.patientunitstayid;

CREATE INDEX ON eicu_time_zero (patientunitstayid);

-- PHASE 2: WBC long format --------------------------------------------------
DROP TABLE IF EXISTS eicu_wbc_long;
CREATE TEMP TABLE eicu_wbc_long AS
SELECT
    s.patientunitstayid,
    l.labresultoffset AS minutes_from_admit,
    (l.labresultoffset - tz.time_zero_min) / 60.0 AS hrs,
    l.labresult AS wbc
FROM eicu_icu.lab l
JOIN eicu_sepsis_cohort s ON l.patientunitstayid = s.patientunitstayid
JOIN eicu_time_zero tz ON s.patientunitstayid = tz.patientunitstayid
WHERE l.labname = 'WBC x 1000'
  AND l.labresult IS NOT NULL
  AND l.labresult > 0
  AND l.labresult < 100
  AND (l.labresultoffset - tz.time_zero_min) / 60.0 BETWEEN -8 AND 96;

CREATE INDEX ON eicu_wbc_long (patientunitstayid);

-- PHASE 3: Lactate long format ----------------------------------------------
DROP TABLE IF EXISTS eicu_lactate_long;
CREATE TEMP TABLE eicu_lactate_long AS
SELECT
    s.patientunitstayid,
    l.labresultoffset AS minutes_from_admit,
    (l.labresultoffset - tz.time_zero_min) / 60.0 AS hrs,
    l.labresult AS lactate
FROM eicu_icu.lab l
JOIN eicu_sepsis_cohort s ON l.patientunitstayid = s.patientunitstayid
JOIN eicu_time_zero tz ON s.patientunitstayid = tz.patientunitstayid
WHERE LOWER(l.labname) = 'lactate'
  AND l.labresult IS NOT NULL AND l.labresult > 0 AND l.labresult < 30
  AND (l.labresultoffset - tz.time_zero_min) / 60.0 BETWEEN -8 AND 96;

CREATE INDEX ON eicu_lactate_long (patientunitstayid);

-- PHASE 4: Microbiology -----------------------------------------------------
-- 4a: Culture results
DROP TABLE IF EXISTS eicu_micro_culture;
CREATE TEMP TABLE eicu_micro_culture AS
SELECT DISTINCT
    s.patientunitstayid,
    CASE WHEN COUNT(m.microlabid) > 0 THEN TRUE ELSE FALSE END AS has_culture_data,
    CASE WHEN COUNT(m.microlabid) FILTER (
        WHERE m.organism IS NOT NULL AND m.organism != ''
          AND LOWER(m.organism) !~* '(no growth|no organism|sterile|normal flora|contaminant)'
    ) > 0 THEN TRUE ELSE FALSE END AS culture_positive,
    BOOL_OR(CASE WHEN LOWER(m.organism) ~* '(escherichia|klebsiella|enterobacter|serratia|pseudomonas|acinetobacter|proteus|haemophilus|neisseria|bacteroides)'
                 THEN TRUE ELSE FALSE END) AS cult_gn,
    BOOL_OR(CASE WHEN LOWER(m.organism) ~* '(staphylococcus|streptococcus|enterococcus|listeria|clostridium)'
                 THEN TRUE ELSE FALSE END) AS cult_gp,
    BOOL_OR(CASE WHEN LOWER(m.organism) ~* '(candida|aspergillus|cryptococcus|yeast|fungus)'
                 THEN TRUE ELSE FALSE END) AS cult_fungal
FROM eicu_sepsis_cohort s
LEFT JOIN eicu_icu.microlab m ON s.patientunitstayid = m.patientunitstayid
GROUP BY s.patientunitstayid;

-- 4b: Clinical diagnosis with pathogen coding
DROP TABLE IF EXISTS eicu_micro_dx;
CREATE TEMP TABLE eicu_micro_dx AS
SELECT DISTINCT
    s.patientunitstayid,
    TRUE AS has_pathogen_dx,
    BOOL_OR(CASE WHEN LOWER(d.diagnosisstring) ~* 'gram.negative|escherichia|klebsiella|pseudomonas|enterobacter|proteus|acinetobacter|haemophilus|neisseria'
                 THEN TRUE ELSE FALSE END) AS dx_gn,
    BOOL_OR(CASE WHEN LOWER(d.diagnosisstring) ~* 'gram.postive|staphylococcus|streptococcus|enterococcus|listeria'
                 THEN TRUE ELSE FALSE END) AS dx_gp,
    BOOL_OR(CASE WHEN LOWER(d.diagnosisstring) ~* 'candida|fungal|aspergillus'
                 THEN TRUE ELSE FALSE END) AS dx_fungal,
    BOOL_OR(CASE WHEN LOWER(d.diagnosisstring) ~* 'bacteremia'
                 THEN TRUE ELSE FALSE END) AS dx_bacteremia
FROM eicu_sepsis_cohort s
JOIN eicu_icu.diagnosis d ON s.patientunitstayid = d.patientunitstayid
WHERE LOWER(d.diagnosisstring) ~* '(bacteremia|gram|staphylococcus|streptococcus|enterococcus|candida|klebsiella|pseudomonas|escherichia|fungal)'
GROUP BY s.patientunitstayid;

-- 4c: Combined
DROP TABLE IF EXISTS eicu_micro;
CREATE TEMP TABLE eicu_micro AS
SELECT
    c.patientunitstayid,
    COALESCE(mc.has_culture_data, FALSE) AS has_culture,
    COALESCE(mc.culture_positive, FALSE) AS culture_pos,
    COALESCE(mc.cult_gn, FALSE) OR COALESCE(md.dx_gn, FALSE) AS g_negative,
    COALESCE(mc.cult_gp, FALSE) OR COALESCE(md.dx_gp, FALSE) AS g_positive,
    COALESCE(mc.cult_fungal, FALSE) OR COALESCE(md.dx_fungal, FALSE) AS fungal,
    COALESCE(md.has_pathogen_dx, FALSE) AS has_pathogen_dx
FROM eicu_sepsis_cohort c
LEFT JOIN eicu_micro_culture mc ON c.patientunitstayid = mc.patientunitstayid
LEFT JOIN eicu_micro_dx md ON c.patientunitstayid = md.patientunitstayid;

-- PHASE 5: Antibiotic spectrum ----------------------------------------------
DROP TABLE IF EXISTS eicu_abx_spectrum;
CREATE TEMP TABLE eicu_abx_spectrum AS
SELECT
    c.patientunitstayid,
    BOOL_OR(LOWER(m.drugname) ~* '(vancomycin|teicoplanin|daptomycin|linezolid)') AS has_gp_abx,
    BOOL_OR(LOWER(m.drugname) ~* '(meropenem|imipenem|ertapenem|cefepime|ceftazidime|piperacillin|aztreonam|ciprofloxacin|levofloxacin|amikacin|gentamicin|tobramycin|colistin|polymyxin)') AS has_gn_abx,
    BOOL_OR(LOWER(m.drugname) ~* '(caspofungin|micafungin|anidulafungin|fluconazole|voriconazole|amphotericin)') AS has_fungal_abx,
    COUNT(DISTINCT m.drugname) AS n_unique_abx
FROM eicu_sepsis_cohort c
JOIN eicu_icu.medication m ON c.patientunitstayid = m.patientunitstayid
WHERE m.drugstartoffset BETWEEN 0 AND 5760  -- first 4 days (96h)
GROUP BY c.patientunitstayid;

-- PHASE 6: Export -----------------------------------------------------------




-- Summary
SELECT
    'eICU Sepsis Cohort' AS label,
    COUNT(*) AS n_stays,
    ROUND(AVG(age)::numeric, 1) AS mean_age,
    ROUND(AVG(CASE WHEN hosp_death = 1 THEN 1.0 ELSE 0.0 END)*100, 1) AS pct_mortality,
    SUM(CASE WHEN has_culture THEN 1 ELSE 0 END) AS n_with_culture,
    SUM(CASE WHEN culture_pos THEN 1 ELSE 0 END) AS n_culture_pos,
    SUM(CASE WHEN has_pathogen_dx THEN 1 ELSE 0 END) AS n_with_pathogen_dx,
    SUM(CASE WHEN g_negative THEN 1 ELSE 0 END) AS n_gn,
    SUM(CASE WHEN g_positive THEN 1 ELSE 0 END) AS n_gp,
    SUM(CASE WHEN fungal THEN 1 ELSE 0 END) AS n_fungal,
    ROUND(AVG(CASE WHEN g_negative OR g_positive OR fungal THEN 1.0 ELSE 0.0 END)*100, 1) AS pct_any_pathogen
FROM (
    SELECT c.*, m.*
    FROM eicu_sepsis_cohort c
    LEFT JOIN eicu_micro m ON c.patientunitstayid = m.patientunitstayid
) sub;

-- ============================================================================
-- EXPORT
\copy (SELECT c.patientunitstayid, c.age, c.gender, c.unittype, c.hosp_death, tz.time_zero_min, COALESCE(m.has_culture, FALSE) AS has_culture, COALESCE(m.culture_pos, FALSE) AS culture_pos, COALESCE(m.g_negative, FALSE) AS g_negative, COALESCE(m.g_positive, FALSE) AS g_positive, COALESCE(m.fungal, FALSE) AS fungal, COALESCE(m.has_pathogen_dx, FALSE) AS has_pathogen_dx, COALESCE(abx.has_gp_abx, FALSE) AS gp_abx, COALESCE(abx.has_gn_abx, FALSE) AS gn_abx, COALESCE(abx.has_fungal_abx, FALSE) AS fungal_abx, COALESCE(abx.n_unique_abx, 0) AS n_abx, COALESCE(lac.max_lactate, NULL) AS max_lactate, COALESCE(lac.lactate_measurements, 0) AS n_lactate FROM eicu_sepsis_cohort c LEFT JOIN eicu_time_zero tz ON c.patientunitstayid = tz.patientunitstayid LEFT JOIN eicu_micro m ON c.patientunitstayid = m.patientunitstayid LEFT JOIN eicu_abx_spectrum abx ON c.patientunitstayid = abx.patientunitstayid LEFT JOIN (SELECT patientunitstayid, MAX(lactate) AS max_lactate, COUNT(*) AS lactate_measurements FROM eicu_lactate_long GROUP BY patientunitstayid) lac ON c.patientunitstayid = lac.patientunitstayid) TO 'data/eicu_cohort.csv' CSV HEADER
\copy (SELECT patientunitstayid, hrs, wbc FROM eicu_wbc_long ORDER BY patientunitstayid, hrs, minutes_from_admit, wbc) TO 'data/eicu_blood_long.csv' CSV HEADER
\copy (SELECT patientunitstayid, hrs, lactate FROM eicu_lactate_long ORDER BY patientunitstayid, hrs, minutes_from_admit, lactate) TO 'data/eicu_lactate_long.csv' CSV HEADER
\copy (SELECT patientunitstayid, organism, culturesite, culturetakenoffset FROM eicu_icu.microlab WHERE patientunitstayid IN (SELECT patientunitstayid FROM eicu_sepsis_cohort) AND organism IS NOT NULL AND organism != '' ORDER BY patientunitstayid, culturetakenoffset, organism, culturesite) TO 'data/eicu_microlab_raw.csv' WITH CSV HEADER

-- ============================================================================
-- eICU STROBE flow: progressive cohort counts for Figure S1
-- ============================================================================

-- ── Step 1: all ICU stays ───────────────────────────────────────────────────
DROP TABLE IF EXISTS es_s1;
CREATE TEMP TABLE es_s1 AS SELECT COUNT(*)::int AS n FROM eicu_icu.patient;

-- ── Step 2: age >= 18 (eICU codes '> 89' and '' specially) ─────────────────
DROP TABLE IF EXISTS es_age;
CREATE TEMP TABLE es_age AS
SELECT p.patientunitstayid, p.unitdischargeoffset, p.age
FROM eicu_icu.patient p
WHERE CASE WHEN p.age = '> 89' THEN 90
           WHEN p.age = '' OR p.age IS NULL THEN NULL
           ELSE CAST(p.age AS INTEGER) END >= 18;
CREATE INDEX ON es_age (patientunitstayid);
DROP TABLE IF EXISTS es_s2;
CREATE TEMP TABLE es_s2 AS SELECT COUNT(*)::int AS n FROM es_age;

-- ── Step 3: ICU length of stay >= 24 h ─────────────────────────────────────
DROP TABLE IF EXISTS es_los;
CREATE TEMP TABLE es_los AS
SELECT a.patientunitstayid
FROM es_age a
WHERE a.unitdischargeoffset >= 1440;
CREATE INDEX ON es_los (patientunitstayid);
DROP TABLE IF EXISTS es_s3;
CREATE TEMP TABLE es_s3 AS SELECT COUNT(*)::int AS n FROM es_los;

-- ── Step 4: sepsis diagnosis (keyword / ICD-9 match) ───────────────────────
DROP TABLE IF EXISTS es_dx;
CREATE TEMP TABLE es_dx AS
SELECT DISTINCT patientunitstayid FROM eicu_icu.diagnosis
WHERE LOWER(diagnosisstring) LIKE '%sepsis%'
   OR LOWER(diagnosisstring) LIKE '%septic%'
   OR LOWER(diagnosisstring) LIKE '%bacteremia%'
   OR icd9code LIKE '038%'
   OR icd9code LIKE '995.91%'
   OR icd9code LIKE '995.92%';
CREATE INDEX ON es_dx (patientunitstayid);

DROP TABLE IF EXISTS es_sepsis;
CREATE TEMP TABLE es_sepsis AS
SELECT l.patientunitstayid
FROM es_los l
JOIN es_dx d ON d.patientunitstayid = l.patientunitstayid;
CREATE INDEX ON es_sepsis (patientunitstayid);
DROP TABLE IF EXISTS es_s4;
CREATE TEMP TABLE es_s4 AS SELECT COUNT(*)::int AS n FROM es_sepsis;

-- ── Step 5: IV antibiotic time zero inside -480..+1440 min ─────────────────
DROP TABLE IF EXISTS es_tz;
CREATE TEMP TABLE es_tz AS
SELECT m.patientunitstayid, MIN(m.drugstartoffset) AS tz_min
FROM eicu_icu.medication m
WHERE LOWER(m.drugname) ~* '(vancomycin|piperacillin|ceftriaxone|cefepime|meropenem|ampicillin|aztreonam|ciprofloxacin|levofloxacin|moxifloxacin|metronidazole|clindamycin|linezolid|ertapenem|imipenem|ceftazidime|cefazolin|gentamicin|tobramycin|amikacin)'
  AND LOWER(m.routeadmin) IN ('iv','iv piggyback','intravenous','iv push')
  AND m.drugstartoffset >= -480
  AND m.drugstartoffset <= 1440
GROUP BY m.patientunitstayid;
CREATE INDEX ON es_tz (patientunitstayid);

DROP TABLE IF EXISTS es_tzok;
CREATE TEMP TABLE es_tzok AS
SELECT s.patientunitstayid, t.tz_min
FROM es_sepsis s
JOIN es_tz t ON t.patientunitstayid = s.patientunitstayid;
CREATE INDEX ON es_tzok (patientunitstayid);
DROP TABLE IF EXISTS es_s5;
CREATE TEMP TABLE es_s5 AS SELECT COUNT(*)::int AS n FROM es_tzok;

-- ── Step 6 (FINAL): >=1 WBC in the -8..+96 h window around time zero ───────
DROP TABLE IF EXISTS es_wbc;
CREATE TEMP TABLE es_wbc AS
SELECT l.patientunitstayid, l.labresultoffset, l.labresult
FROM eicu_icu.lab l
WHERE l.labname = 'WBC x 1000'
  AND l.labresult IS NOT NULL
  AND l.labresult > 0
  AND l.labresult < 100;
CREATE INDEX ON es_wbc (patientunitstayid);

DROP TABLE IF EXISTS es_s6;
CREATE TEMP TABLE es_s6 AS
SELECT COUNT(DISTINCT t.patientunitstayid)::int AS n
FROM es_tzok t
JOIN es_wbc w ON w.patientunitstayid = t.patientunitstayid
WHERE (w.labresultoffset - t.tz_min) / 60.0 BETWEEN -8 AND 96;

-- -- Funnel (written to disk for the R side to read) ----------------------------
\copy (SELECT (SELECT n FROM es_s1) AS step1_icu_stays, (SELECT n FROM es_s2) AS step2_age18, (SELECT n FROM es_s3) AS step3_los24h, (SELECT n FROM es_s4) AS step4_sepsis, (SELECT n FROM es_s5) AS step5_time_zero, (SELECT n FROM es_s6) AS step6_wbc_final) TO 'data/eicu_strobe_counts.csv' CSV HEADER
