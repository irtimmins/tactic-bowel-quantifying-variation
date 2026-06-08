# Synthetic colon-cancer waiting-times cohort

R code to generate a **synthetic** dataset that mirrors the structure of the real
ICON `colon_cohort_2015_2022.rds` (the 172-column cohort produced by the
`Colon_cancer_prepare_data1a.R` + `1b.R` pipeline: NCRAS → HES-APC → CWT linkage).

The synthetic `.dta` lets you develop and debug the Stata waiting-times analysis
**off** the secure server, then run the *same* code unchanged on the real ICON
data. No real patient data is used or required to build the synthetic file.

---

## Two-stage, disclosure-safe workflow

```
  ON THE ICON SERVER                         ANYWHERE (laptop, off-server)
  ------------------                         -----------------------------
  real colon_cohort_2015_2022.rds
            |
            v
  01_extract_icon_metadata.R
   - validates the 172-col structure
   - writes column_spec.csv         --->  (disclosure check)
   - writes icon_metadata.rds  ----------------> 02_build_synthetic_cohort.R
     (AGGREGATE numbers only:                     - simulates 172 cols
      proportions, means, variances)              - enforces consistency rules
   - writes icon_metadata_summary.txt             - writes synthetic_colon_cohort.dta
            |                                                |
            v                                                v
   (small cells suppressed,                        develop & test Stata code,
    counts rounded — safe to export)               then run it on the real .rds/.dta
```

The only thing that crosses the boundary is `icon_metadata.rds` — a list of
**aggregate** summary statistics (category proportions, waiting-time means/variances,
between-hospital variance, year/route/stage effects). Cells smaller than
`MIN_CELL` (default 10) are suppressed and all counts are rounded, so it is
designed to pass a standard ONS/NHS disclosure review. No row-level data is ever
written by Stage 1.

You do **not** have to run Stage 1 first: `02_build_synthetic_cohort.R` runs
standalone using built-in default parameters (taken from the published Stata
logs / provider summaries), so you can produce a usable `.dta` today. Run Stage 1
later to make the synthetic distributions match the real cohort more closely.

---

## Files

| File | Where it runs | What it does |
|------|---------------|--------------|
| `01_extract_icon_metadata.R` | ICON server | Validates structure; writes `column_spec.csv`, `icon_metadata.rds`, `icon_metadata_summary.txt`. |
| `02_build_synthetic_cohort.R` | anywhere | Generates `synthetic_colon_cohort.dta` (+ `.rds`) with the exact 172-column structure. Uses `icon_metadata.rds` if supplied, else defaults. |
| `synthetic_colon_cohort.dta` / `.rds` | — | A ready-to-use default build (50,000 patients, 148 hospitals) so you can start the Stata code immediately. |
| `example_column_spec.csv` | — | Example structural dictionary (generated from the synthetic build — illustrates the format Stage 1 produces). |
| `example_icon_metadata_summary.txt` | — | Example of the human-readable parameter summary. |

### How to run

Stage 1 (on the server), after editing the paths at the top:
```r
source("01_extract_icon_metadata.R")
```
Review `column_spec.csv` and `icon_metadata_summary.txt`, clear disclosure, copy
`icon_metadata.rds` out.

Stage 2 (anywhere):
```r
# standalone (defaults):
source("02_build_synthetic_cohort.R")

# or data-informed: set META_PATH <- "icon_metadata.rds" at the top first
```

Dependencies: **base R only**. `.dta` is written with `haven` if installed,
otherwise with `foreign::write.dta` (note: the `foreign` fallback abbreviates
long variable names, so install `haven` for a faithful `.dta`).

---

## What is modelled, and how faithfully

The 172 columns fall into two groups.

**Analytic core** (≈30 columns) — modelled with realistic distributions and the
relationships the analysis depends on:

- Outcomes `wt_dx_to_dtt`, `wt_dtt_to_tx`, `wt_dx_to_tx` (diagnosis→DTT, DTT→surgery,
  diagnosis→surgery), drawn from skewed (Gamma) distributions with a **hospital
  random effect** so there is genuine between-hospital variation (a non-trivial ICC)
  for the multilevel models, plus modest covariate effects (route, stage, CCI, age,
  calendar year) so the Stata models recover sensible signs.
- Patient/tumour/pathway factors: `agediag`, `sex`, `ethnicity_group_broad`,
  `NHSE_reversed_imd_quintile_lsoas`, `stage`/`stage_best`, `route_combined`,
  `grade`, `colon_proc_type_primary`, `modality`, `ydiag`/`diagmdy`.
- Hospital/trust structure: `diag_hosp` (5-char site) nested within `diag_trust`
  (3-char), `SITETRET`/`PROCODE3`/`org_treat_start` for trust-of-treatment, with a
  minority of patients diagnosed and treated at different trusts.

**Structural filler** (the remaining CWT/HES administrative columns) — generated
with the **correct type and plausible category levels** but simplified
distributions (often a single dominant value or `NA`). They exist so the dataset
has the exact shape and so any code that touches them runs; they are not intended
to reproduce real administrative distributions. If you later need realistic
distributions for any of these, Stage 1 already captures generic summaries for
every column in `column_spec.csv`, and the generator is easy to extend.

### Built-in consistency rules (so derived/QA variables behave correctly)

- `wt_dx_to_tx == wt_dx_to_dtt + wt_dtt_to_tx` wherever DTT is observed.
- Dates regenerated from the simulated intervals: `dtt_date = diagmdy + wt_dx_to_dtt`,
  `tx_date = diagmdy + wt_dx_to_tx`, `EPISTART = ADMIDATE = tx_date`; `mdt_date`
  before `dtt_date`; `cwt_tx_date` within ±5 days of `tx_date` (`tx_date_diff ≤ 5`).
- Pathway ordering Dx ≤ MDT ≤ DTT ≤ Tx, so the `dx_le_*`, `dtt_le_tx_ok`, `seq_ok`
  flags compute to the same logic as `1b`.
- `days_diag_to_surg = wt_dx_to_tx`, bounded to **[0, 180]** exactly as the real
  cohort's `1b` filter enforces.
- `agediag`/`birthmdy`/`diagmdy`/`ydiag` mutually consistent; `m_best = "M0"`,
  `stage ∈ {1,2,3}`, `route_combined` excludes Emergency presentation / Unknown,
  `modality ∈ {01,23,24}` — matching the cohort's inclusion criteria.
- ~10% of records have a **missing DTT** (`dtt_date`, `wt_dx_to_dtt`,
  `wt_dtt_to_tx` all `NA`) while keeping `wt_dx_to_tx` — the real cohort only
  requires `wt_dx_to_tx ≥ 0`, so this reproduces that pattern and exercises the
  Stata `!missing(...)` filters.

---

## Notes for the Stata side (next step)

- **Variable-name bridge.** Your patient-level analysis (`Analysis_tactic_rapid.do`)
  uses the *rapid* names (`age_rapid`, `sex_rapid`, `stage_rapid`, `time_diag_dtt`, …)
  while the provider-level/random-effects code uses the *ICON* names
  (`agediag`, `sex`, `stage`, `wt_dx_to_dtt`, …). This synthetic file uses the
  **ICON names** (it mirrors the ICON `.rds`). Plan a short `rename`/harmonisation
  block at the top of the Stata code so the same models run on both the synthetic
  file, the real ICON `.dta`, and your colleague's rapid data.
- **`laterality.x` / `.y`.** Stata variable names cannot contain a dot, so the
  `.dta` writer converts these to `laterality_x` / `laterality_y` (the same
  conversion Stata applies on import). The `.rds` keeps the exact original names.
- **Provider-level covariates** (`comprehensive`, `teaching`, `cqc_rating`,
  `staff_eng_cat`, `moral_cat`, `bed_occ_cat`, `competitor_status`) and
  **travel times** are *not* part of the 172-column cohort — they are merged in
  later from external files (`NHSHospitals_services_*.xlsx`, the drive-time matrix,
  the net-gain file). If you want to test that merge too, the next extension is a
  small generator that emits one synthetic row per `diag_hosp` with those columns;
  flag it and it can be added.

## Tunable parameters (top of `02_build_synthetic_cohort.R`)

`N_PATIENTS`, `N_HOSPITALS`, `SEED`, `OUT_DIR`, `META_PATH`, plus the `defaults`
list (year weights, waiting-time means/shapes, between-hospital SDs, category
probabilities, covariate effect sizes, missingness rates).
