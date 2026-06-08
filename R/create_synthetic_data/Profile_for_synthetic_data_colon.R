# =============================================================================
# Colon waiting times - profiling for synthetic data generation
# -----------------------------------------------------------------------------
# Runs on the secure server. Extracts the aggregate distributions needed to
# generate a synthetic NCRAS+HES cohort and a synthetic CWT records table that
# reproduce the cancer-waiting-times merge (the 1b linkage step) and its
# patient-exclusion funnel.
#
# Disclosure control: counts rounded to nearest 5, cells < 10 suppressed,
# intervals reported as quantiles only. Output is aggregate. Still send the
# saved objects through output checking before they leave the environment.
#
# Produces: colon_profile_for_synthetic.rds   (distributions)
#           colon_pipeline_spec.rds            (column spec + merge constants)
# Inputs:   colon_cohort_2015_2022.rds         (post-merge analysis cohort)
#           colon_cohort_cci_2015_2022.rds      (optional, for cci_group)
#           the partitioned CWT dataset         (optional, per-record stats)
# =============================================================================

library(tidyverse)
library(arrow)
library(lubridate)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"
save_dir <- "Data/synthetic/"
cwt_path <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/11_CWT_data_partitioned"
SDC_MIN  <- 10L
colon_icd <- c("C18","C180","C181","C182","C183","C184","C185","C186","C187","C188","C189")

cohort <- readRDS(paste0(base_dir, "colon_cohort_cci_2015_2022.rds"))

# -----------------------------------------------------------------------------
# Normalise cohort variables to canonical forms before profiling.
# The real cohort may store these with labelled factors or abbreviations that
# would flow through to the profile and cause mismatches in the generator.
# -----------------------------------------------------------------------------
cohort <- cohort %>%
  mutate(
    # imd: labelled factor levels like "1 - most deprived" -> integer 1..5
    NHSE_reversed_imd_quintile_lsoas = as.integer(substr(
      as.character(NHSE_reversed_imd_quintile_lsoas), 1, 1
    )),
    # screendetected: Y/N/blank -> 1/0
    screendetected = case_when(
      as.character(screendetected) == "Y" ~ 1L,
      as.character(screendetected) == "N" ~ 0L,
      TRUE                               ~ 0L
    ),
    # dco: Y/N -> 1/0
    dco = case_when(
      as.character(dco) == "Y" ~ 1L,
      as.character(dco) == "N" ~ 0L,
      TRUE                     ~ 0L
    ),
    # route_combined: normalise any abbreviated label to the full form
    route_combined = recode(as.character(route_combined), "TWW" = "Two Week Wait")
  )

# -----------------------------------------------------------------------------
# Disclosure-safe helpers
# -----------------------------------------------------------------------------
r5  <- function(x) round(x / 5) * 5
sup <- function(n) ifelse(n < SDC_MIN, NA_real_, r5(n))

cat_marg <- function(df, var) {
  if (!var %in% names(df)) return(NULL)
  df %>%
    count(.data[[var]], name = "n") %>%
    mutate(n_safe = sup(n), prop = round(n / sum(n), 4)) %>%
    select(level = 1, n_safe, prop)
}

q_sum <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < SDC_MIN) return(NULL)
  tibble(
    n_nonmiss = r5(length(x)),
    p05 = quantile(x, .05), p10 = quantile(x, .10), p25 = quantile(x, .25),
    p50 = quantile(x, .50), p75 = quantile(x, .75), p90 = quantile(x, .90),
    p95 = quantile(x, .95),
    mean = round(mean(x), 1), sd = round(sd(x), 1)
  )
}

profile <- list()

# =============================================================================
# A. Scale and temporal structure
# =============================================================================
profile$n_patients <- r5(nrow(cohort))
profile$by_year    <- cat_marg(cohort, "ydiag")

# =============================================================================
# B. Patient, tumour and pathway marginals (drive covariate sampling)
# =============================================================================
cohort2 <- cohort %>%
  mutate(age_grp = cut(agediag,
                       breaks = c(-Inf,50,55,60,65,70,75,80,85,Inf),
                       labels = c("<50","50-54","55-59","60-64","65-69",
                                  "70-74","75-79","80-84","85+"),
                       right = FALSE))

profile$marginals <- list(
  sex             = cat_marg(cohort,  "sex"),
  age_grp         = cat_marg(cohort2, "age_grp"),
  ethnicity       = cat_marg(cohort,  "ethnicity_group_broad"),
  imd_quintile    = cat_marg(cohort,  "NHSE_reversed_imd_quintile_lsoas"),
  stage           = cat_marg(cohort,  "stage"),
  route_combined  = cat_marg(cohort,  "route_combined"),
  screendetected  = cat_marg(cohort,  "screendetected"),
  proc_type       = cat_marg(cohort,  "colon_proc_type_primary"),
  dco             = cat_marg(cohort,  "dco"),
  emergency       = cat_marg(cohort,  "emergency"),
  cci_group       = cat_marg(cohort,  "cci_group"),       # NULL if absent
  perf_status     = cat_marg(cohort,  "perf_status")      # NULL if absent
)

# stage x year (captures any COVID-era stage shift)
profile$stage_by_year <- cohort %>%
  count(ydiag, stage, name = "n") %>%
  group_by(ydiag) %>%
  mutate(prop = round(n / sum(n), 4), n_safe = sup(n)) %>%
  ungroup() %>%
  select(ydiag, stage, n_safe, prop)

# =============================================================================
# C. Waiting-time intervals: overall, and the between-hospital signal
# =============================================================================
profile$intervals_overall <- list(
  wt_dx_to_dtt = q_sum(cohort$wt_dx_to_dtt[cohort$wt_dx_to_dtt >= 0]),
  wt_dtt_to_tx = q_sum(cohort$wt_dtt_to_tx[cohort$wt_dtt_to_tx >= 0]),
  wt_dx_to_tx  = q_sum(cohort$wt_dx_to_tx[cohort$wt_dx_to_tx  >= 0]),
  surv_from_dx_days = q_sum(as.integer(cohort$finmdy - cohort$diagmdy))
)

# between-hospital SD of the per-hospital mean wait (random-intercept signal)
hosp_sd <- function(num, hosp) {
  d <- tibble(num = num, hosp = hosp) %>%
    filter(!is.na(num), num >= 0, grepl("^R[A-Z0-9]{4}$", hosp)) %>%
    group_by(hosp) %>% filter(n() >= SDC_MIN) %>%
    summarise(m = mean(num), .groups = "drop")
  if (nrow(d) < 2) return(NA_real_)
  round(sd(d$m), 2)
}
profile$between_hosp_sd <- list(
  wt_dx_to_dtt = hosp_sd(cohort$wt_dx_to_dtt, cohort$diag_hosp),
  wt_dtt_to_tx = hosp_sd(cohort$wt_dtt_to_tx, cohort$diag_hosp),
  wt_dx_to_tx  = hosp_sd(cohort$wt_dx_to_tx,  cohort$diag_hosp)
)

# =============================================================================
# D. Trust / hospital volume structure and change-of-trust rate
# =============================================================================
vol <- function(v) {
  v <- v[!is.na(v) & v != ""]
  if (!length(v)) return(NULL)
  sizes <- as.integer(table(v))
  tibble(n_distinct = length(sizes),
         vol_p25 = quantile(sizes, .25), vol_p50 = quantile(sizes, .50),
         vol_p75 = quantile(sizes, .75), vol_p90 = quantile(sizes, .90))
}
hosp_qc <- cohort$diag_hosp[grepl("^R[A-Z0-9]{4}$", cohort$diag_hosp)]
profile$volume <- list(
  diag_trust = vol(cohort$diag_trust),
  diag_hosp  = vol(hosp_qc),
  PROCODE3   = vol(cohort$PROCODE3)
)

# diagnosis trust vs treatment trust (3-char), if both present
if (all(c("diag_trust","PROCODE3") %in% names(cohort))) {
  profile$change_trust_rate <- cohort %>%
    filter(!is.na(diag_trust), !is.na(PROCODE3)) %>%
    summarise(rate = round(mean(substr(diag_trust,1,3) != substr(PROCODE3,1,3)), 4)) %>%
    pull(rate)
}

# =============================================================================
# E. CWT per-record structure from the raw partitioned dataset
# =============================================================================
if (dir.exists(cwt_path)) {
  ids <- cohort$pseudo_patientid
  cwt <- open_dataset(cwt_path) %>%
    filter(site_icd10 %in% colon_icd) %>%
    collect() %>%
    mutate(pseudo_patientid = as.character(pseudo_patientid)) %>%
    filter(pseudo_patientid %in% ids)

  profile$cwt_records_per_patient <- cwt %>%
    count(pseudo_patientid, name = "k") %>%
    count(k, name = "n_pat") %>%
    mutate(n_safe = sup(n_pat), prop = round(n_pat / sum(n_pat), 4)) %>%
    select(records = k, n_safe, prop)

  profile$cwt_modality  <- cat_marg(cwt, "modality")
  profile$cwt_site_icd10 <- cat_marg(cwt, "site_icd10")

  profile$cwt_coverage <- tibble(
    pct_any_cwt = round(mean(ids %in% cwt$pseudo_patientid), 4)
  )
} else {
  profile$cwt_records_per_patient <- NULL
  profile$cwt_modality  <- NULL
  profile$cwt_site_icd10 <- NULL
  profile$cwt_coverage  <- NULL
}

# mdt completeness and timing relative to dtt, from the post-merge cohort
profile$cwt_completeness <- tibble(
  pct_mdt = round(mean(!is.na(cohort$mdt_date)), 4),
  pct_dtt = round(mean(!is.na(cohort$dtt_date)), 4)
)
profile$mdt_to_dtt <- cohort %>%
  filter(!is.na(mdt_date), !is.na(dtt_date)) %>%
  mutate(d = as.integer(dtt_date - mdt_date)) %>% pull(d) %>% q_sum()

# =============================================================================
# F. Merge-glue distributions (signed) - reproduce the linkage behaviour
#    first treatment = HES surgery date (tx_date); cwt_tx_date = treat_start
# =============================================================================
glue <- cohort %>%
  filter(!is.na(dtt_date)) %>%
  mutate(
    days_dx_to_dtt   = as.integer(dtt_date - diagmdy),
    dtt_to_cwt_treat = as.integer(cwt_tx_date - dtt_date),
    cwt_vs_tx        = as.integer(cwt_tx_date - tx_date)
  )
profile$cwt_glue <- list(
  days_dx_to_dtt   = q_sum(glue$days_dx_to_dtt),
  dtt_to_cwt_treat = q_sum(glue$dtt_to_cwt_treat),
  cwt_vs_tx        = q_sum(glue$cwt_vs_tx)
)
profile$cwt_agreement <- cohort %>%
  filter(!is.na(cwt_tx_date)) %>%
  mutate(cwt_vs_tx = as.integer(cwt_tx_date - tx_date)) %>%
  summarise(
    pct_exact    = round(mean(cwt_vs_tx == 0), 4),
    pct_within_5 = round(mean(abs(cwt_vs_tx) <= 5), 4)
  )

# share of linked patients with a missing DTT (kept on wt_dx_to_tx alone)
profile$pct_dtt_missing_in_cohort <- round(mean(is.na(cohort$dtt_date)), 4)

# tx_date_diff distribution on the final cohort (all <= 5 by construction)
profile$tx_date_diff <- q_sum(cohort$tx_date_diff)

# =============================================================================
# Pipeline spec: column manifest + merge constants for the generator/validator
# =============================================================================
spec_type <- function(x) {
  if (inherits(x, "Date")) "Date" else if (is.factor(x)) "factor" else
    if (is.logical(x)) "logical" else if (is.integer(x)) "integer" else
      if (is.numeric(x)) "numeric" else "character"
}

cohort_spec <- tribble(
  ~name,                              ~tier,
  "pseudo_patientid",                 "required",
  "pseudo_tumourid",                  "core",
  "diagmdy",                          "required",
  "ydiag",                            "required",
  "yearmonth_diag",                   "core",
  "sitestr",                          "required",
  "stage",                            "required",
  "stage_best",                       "core",
  "agediag",                          "required",
  "sex",                              "required",
  "ethnicity_group_broad",            "required",
  "NHSE_reversed_imd_quintile_lsoas", "required",
  "route_combined",                   "required",
  "screendetected",                   "core",
  "cci_group",                        "core",
  "rcs_ch_score",                     "required",
  "perf_status",                      "core",
  "colon_proc_type_primary",          "core",
  "any_mets",                         "core",
  "tnm_m",                            "core",
  "pretreat_m",                       "core",
  "diag_hosp",                        "required",
  "diag_trust",                       "required",
  "first_trust",                      "core",
  "tx_date",                          "required",
  "days_diag_to_surg",                "required",
  "tx_hosp",                          "required",
  "tx_trust",                         "required",
  "SITETRET",                         "core",
  "PROCODE3",                         "core",
  "emergency",                        "core",
  "dco",                              "core",
  "dead",                             "core",
  "finmdy",                           "core",
  "change_trust",                     "core",
  # filled by the CWT merge
  "site_icd10",                       "required",
  "modality",                         "required",
  "dtt_date",                         "required",
  "cwt_tx_date",                      "required",
  "mdt_date",                         "core",
  "tx_date_diff",                     "required",
  "diff_cwt_cr_treat",                "required",
  "diff_cwt_cr_treat_cat",            "required",
  "missing_CWT",                      "required",
  "tx_trust_cwt",                     "required",
  "site_match",                       "required",
  "dx_le_mdt_ok",                     "core",
  "dx_le_dtt_ok",                     "core",
  "dx_le_tx_ok",                      "core",
  "mdt_le_dtt_ok",                    "core",
  "dtt_le_tx_ok",                     "core",
  "mdt_le_tx_ok",                     "core",
  "seq_ok",                           "required",
  "wt_dx_to_dtt",                     "required",
  "wt_dtt_to_tx",                     "required",
  "wt_dx_to_tx",                      "required"
)
# attach the observed type for any column that exists in the real cohort
cohort_spec <- cohort_spec %>%
  mutate(type = map_chr(name, ~ if (.x %in% names(cohort)) spec_type(cohort[[.x]]) else NA_character_))

cwt_spec <- tribble(
  ~name,                ~type,        ~tier,
  "pseudo_patientid",   "character",  "required",
  "site_icd10",         "character",  "required",
  "modality",           "character",  "required",
  "org_treat_start",    "character",  "required",
  "crtp_date",          "character",  "core",
  "date_first_seen",    "character",  "core",
  "mdt_date",           "character",  "core",
  "treat_period_start", "character",  "required",
  "treat_start",        "character",  "required"
)

pipeline_spec <- list(
  cohort_spec   = cohort_spec,
  cwt_spec      = cwt_spec,
  colon_icd     = colon_icd,
  stage_levels  = c("1","2","3"),
  merge_const   = list(
    tx_window_days   = c(-90L, 180L),     # NCRAS-HES window (days_diag_to_surg)
    tx_date_diff_max = 5L,                # HES vs CWT treatment-date agreement
    modality_keep    = c("01","23","24"),
    modality_2324_from = as.Date("2020-06-01")
  )
)

saveRDS(profile,       paste0(save_dir, "colon_profile_for_synthetic.rds"))
saveRDS(pipeline_spec, paste0(save_dir, "colon_pipeline_spec.rds"))
cat("Saved colon_profile_for_synthetic.rds (", length(profile), "sections) and colon_pipeline_spec.rds\n")
str(profile, max.level = 1)
