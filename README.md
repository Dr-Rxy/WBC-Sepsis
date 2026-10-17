# WBC-Sepsis

Analysis code (SQL and R) for the manuscript:

> **Risk stratification of sepsis using longitudinal white blood cell trajectories:
> microbiological differences and multi-database external validation**

Derivation cohort **MIMIC-IV v3.1** (n = 26,369); external validation in **eICU-CRD** and
**NWICU** (pooled N = 53,074). The analysis identifies two WBC-trajectory subgroups
(C1 high-WBC vs C2 normal-WBC).

This repository contains analysis code only. No patient-level data are included.

## Contents

| Path | Content |
|---|---|
| `sql/01_mimic.sql` | MIMIC-IV cohort construction and CSV export |
| `sql/02_eicu.sql` | eICU-CRD cohort construction and CSV export |
| `sql/03_nwicu.sql` | NWICU cohort construction and CSV export |
| `R/00_functions.R` | Shared function library (clustering, Cox, mediation, logging) |
| `R/01_cohort_analyses.R` | Derivation cohort and both external validations (Tables 0–8) |
| `R/02_secondary_outputs.R` | Secondary and supplementary analyses |

## Window specification

The four clustering windows are **half-open and mutually disjoint**, and are named by
their **start** hour. Every measurement falls in exactly one window.

| Feature | Interval (hours from time zero) | Width |
|---|---|---|
| `wbc_baseline` | `[-8, 24)` | 32 h |
| `wbc_t24h` | `[24, 48)` | 24 h |
| `wbc_t48h` | `[48, 72)` | 24 h |
| `wbc_t72h` | `[72, 96)` | 24 h |

The baseline window extends 8 hours backwards so that pre-antibiotic values are included.
The within-window median of each feature is the clustering input. The extraction range is
`[-8, 96)` hours. There is deliberately **no separate 96-hour window**: that period is
already occupied by `wbc_t72h`, so a further window would either overlap it or fall partly
outside the extraction range.

> **Note on earlier revisions of this repository.** The first commit used a different grid —
> baseline `[-8, 24]` followed by `[12, 36)`, `[36, 60)`, `[60, 84)`, of which the first pair
> overlaps — and also defined a separate 84–108 h window. The cohort size, the subgroup
> sizes and every effect estimate differ between the two grids, so the earlier revision does
> **not** reproduce the manuscript. Only the current revision does.

## Time zero and the laboratory window

All three cohorts use a single anchor: **time zero is the time of the first intravenous
antibiotic administration**, with a documented fallback where that time cannot be
determined (MIMIC-IV uses the Sepsis-3 time of suspected infection, eICU the time of ICU
admission). Laboratory maxima, lactate included, are taken over **-24 h to +72 h around
time zero** in every cohort.

> **Note on earlier revisions.** The previous revision took the laboratory maxima over
> windows anchored on ICU admission, and used different widths in different databases
> (-8 to +96 h in eICU, to +120 h in NWICU). That put the laboratory values on a different
> clock from the four clustering windows above. The current revision is the one that
> reproduces the manuscript.

## Requirements

- PostgreSQL with the three databases loaded under the schema names used in the SQL files
- `psql` on `PATH`
- R ≥ 4.6.0 with `data.table`, `cluster`, `mclust`, `survival`, `pROC`, `mice` and `mediation`;
  `Rscript` on `PATH`

MIMIC must be **v3.1**: earlier releases yield a different cohort and the results are not
comparable.

## How to run

Run from the repository root — the SQL scripts export to the relative path `data/`.

```bash
mkdir -p data
export PGPASSWORD=<your password>

psql -h localhost -p 5432 -U postgres -d mimic4_v31 -v ON_ERROR_STOP=1 -f sql/01_mimic.sql
psql -h localhost -p 5432 -U postgres -d eicu      -v ON_ERROR_STOP=1 -f sql/02_eicu.sql
psql -h localhost -p 5432 -U postgres -d nwicu     -v ON_ERROR_STOP=1 -f sql/03_nwicu.sql

Rscript R/01_cohort_analyses.R
Rscript R/02_secondary_outputs.R
```

`R/01` must finish before `R/02`. Result tables are written to `output/tables/`.
Expected cohort sizes: 26,369 (MIMIC-IV) / 21,781 (eICU) / 4,924 (NWICU).

## Data availability

No patient-level data are distributed here. MIMIC-IV, eICU-CRD and NWICU are available to
credentialed users via PhysioNet (https://physionet.org) after the required training and data
use agreement.
