# =============================================================================
# Colon waiting times - synthetic data validation
# -----------------------------------------------------------------------------
# Three checks, run off-server on the synthetic files:
#   1. Conformance - the analysis cohort matches the spec (columns, types,
#      tiers, stage levels, unique IDs).
#   2. Internal consistency - re-derive the waiting times and the pathway
#      ordering flags from the stored dates and confirm they match the stored
#      columns, and that the linkage filters hold (tx_date_diff <= 5, total
#      wait in range).
#   3. Synthetic vs profile - quick comparison of key marginals and the
#      exclusion funnel reproduced from Table A + Table B.
#
# Needs colon_pipeline_spec.rds (aggregate, travels with the bundle).
# =============================================================================

library(tidyverse)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"
save_dir <- "Data/synthetic/"
spec     <- readRDS(paste0(save_dir, "colon_pipeline_spec.rds"))
syn      <- readRDS(paste0(save_dir, "colon_cohort_synthetic.rds"))
A        <- readRDS(paste0(save_dir, "colon_ncras_hes_synthetic.rds"))
cwt      <- readRDS(paste0(save_dir, "colon_cwt_records_synthetic.rds"))
mc       <- spec$merge_const

obs_type <- function(x) {
  if (inherits(x, "Date")) "Date" else if (is.factor(x)) "factor" else
    if (is.logical(x)) "logical" else if (is.integer(x)) "integer" else
      if (is.numeric(x)) "numeric" else "character"
}
compat <- function(e, o) is.na(e) || e == o ||
  (e %in% c("integer","numeric") && o %in% c("integer","numeric"))

# -----------------------------------------------------------------------------
# 1. Conformance
# -----------------------------------------------------------------------------
cs <- spec$cohort_spec
present  <- cs$name %in% names(syn)
miss_req <- cs$name[!present & cs$tier == "required"]
miss_core<- cs$name[!present & cs$tier == "core"]

rows <- cs %>% filter(name %in% names(syn))
obs  <- vapply(rows$name, function(nm) obs_type(syn[[nm]]), character(1))
mism <- rows %>% mutate(observed = obs) %>%
  filter(!mapply(compat, type, observed)) %>% select(name, expected = type, observed)

bad_stage <- setdiff(unique(na.omit(syn$stage)), spec$stage_levels)

cat("== 1. Conformance ==\n")
cat("Rows:", nrow(syn), " Cols:", ncol(syn), "\n")
cat("Missing required cols:", if (length(miss_req)) paste(miss_req, collapse=", ") else "none", "\n")
cat("Missing core cols:    ", if (length(miss_core)) paste(miss_core, collapse=", ") else "none", "\n")
if (nrow(mism)) { cat("Type mismatches:\n"); print(mism) } else cat("Type mismatches:     none\n")
cat("Unexpected stage levels:", if (length(bad_stage)) paste(bad_stage, collapse=", ") else "none", "\n")
cat("Duplicate patient IDs:", sum(duplicated(syn$pseudo_patientid)), "\n\n")

# -----------------------------------------------------------------------------
# 2. Internal consistency
# -----------------------------------------------------------------------------
re <- syn %>% mutate(
  re_wt_dx_to_dtt = as.numeric(dtt_date - diagmdy),
  re_wt_dx_to_tx  = as.numeric(tx_date  - diagmdy),
  re_wt_dtt_to_tx = as.numeric(tx_date  - dtt_date),
  re_seq_ok = (is.na(mdt_date) | diagmdy <= mdt_date) &
    (is.na(dtt_date) | diagmdy <= dtt_date) &
    (diagmdy <= tx_date) &
    (is.na(mdt_date) | is.na(dtt_date) | mdt_date <= dtt_date) &
    (is.na(dtt_date) | dtt_date <= tx_date) &
    (is.na(mdt_date) | mdt_date <= tx_date),
  re_site_match = as.integer(str_sub(site_icd10,1,3) == str_sub(sitestr,1,3))
)
eq <- function(a, b) (a == b) | (is.na(a) & is.na(b))

cat("== 2. Internal consistency ==\n")
cat("wt_dx_to_dtt re-derives:", round(100*mean(eq(re$re_wt_dx_to_dtt, re$wt_dx_to_dtt)), 2), "%\n")
cat("wt_dtt_to_tx re-derives:", round(100*mean(eq(re$re_wt_dtt_to_tx, re$wt_dtt_to_tx)), 2), "%\n")
cat("wt_dx_to_tx  re-derives:", round(100*mean(eq(re$re_wt_dx_to_tx,  re$wt_dx_to_tx)),  2), "%\n")
cat("seq_ok       re-derives:", round(100*mean(eq(re$re_seq_ok, re$seq_ok)), 2), "%\n")
cat("site_match   re-derives:", round(100*mean(eq(re$re_site_match, re$site_match)), 2), "%\n")
id_ok <- with(syn, wt_dx_to_tx == wt_dx_to_dtt + wt_dtt_to_tx)
cat("sum identity (where DTT observed):", round(100*mean(id_ok[!is.na(id_ok)]), 2), "%\n")
cat("tx_date_diff <= ", mc$tx_date_diff_max, ":", all(syn$tx_date_diff <= mc$tx_date_diff_max), "\n")
cat("total wait within [0,180]:", all(syn$wt_dx_to_tx >= 0 & syn$wt_dx_to_tx <= 180), "\n")
cat("diag_hosp nested in diag_trust:",
    all(substr(syn$diag_hosp,1,3) == syn$diag_trust), "\n")
cat("tx_trust == substr(tx_hosp,1,3):",
    all(syn$tx_trust == substr(syn$tx_hosp, 1, 3), na.rm = TRUE), "\n")
if (all(c("tx_trust","tx_trust_cwt") %in% names(syn))) {
  agree <- mean(syn$tx_trust == syn$tx_trust_cwt, na.rm = TRUE)
  cat(sprintf("tx_trust (HES) vs tx_trust_cwt (CWT) agreement: %.1f%%  (expect ~95%% by design)\n",
              100 * agree))
}
# diff_cwt_cr_treat_cat: since the cohort is pre-filtered to tx_date_diff <= 5,
# expect only cats 1 (exact), 2 (within 5d), and occasionally 5 (edge)
if ("diff_cwt_cr_treat_cat" %in% names(syn)) {
  cat("\ndiff_cwt_cr_treat_cat distribution (cats 3/4 absent = correct for pre-filtered cohort):\n")
  print(table(syn$diff_cwt_cr_treat_cat, useNA = "ifany"))
}
# quick check of Table A for the new columns added in the last round
a_new <- c("yearmonth_diag","any_mets","tnm_m","pretreat_m","rcs_ch_score")
cat("\nTable A new columns present:", paste(a_new[a_new %in% names(A)], collapse=", "), "\n")
cat("Table A new columns absent: ",
    if (length(setdiff(a_new, names(A)))) paste(setdiff(a_new, names(A)), collapse=", ") else "none", "\n\n")

# -----------------------------------------------------------------------------
# 3. Synthetic vs profile, and the exclusion funnel from the raw tables
# -----------------------------------------------------------------------------
prof_path <- paste0(save_dir, "colon_profile_for_synthetic.rds")
if (file.exists(prof_path)) {
  prof <- readRDS(prof_path)
  show_marg <- function(var, pkey) {
    if (is.null(prof$marginals[[pkey]])) return(invisible())
    cmp <- syn %>% count(.data[[var]], name = "n_syn") %>%
      mutate(prop_syn = round(n_syn / sum(n_syn), 3), level = as.character(.data[[var]])) %>%
      left_join(prof$marginals[[pkey]] %>% mutate(level = as.character(level)) %>%
                  select(level, prop_real = prop), by = "level") %>%
      select(level, prop_syn, prop_real)
    cat("--", var, "--\n"); print(cmp); cat("\n")
  }
  cat("== 3. Marginals: synthetic vs profile ==\n")
  show_marg("stage", "stage")
  show_marg("route_combined", "route_combined")
  show_marg("sex", "sex")
}

# Reproduce the funnel directly from Table A + Table B (the same 1b logic)
cwt_p <- cwt %>%
  mutate(dtt_date = as.Date(treat_period_start, "%d/%m/%Y"),
         cwt_tx_date = as.Date(treat_start, "%d/%m/%Y")) %>%
  filter(modality %in% mc$modality_keep,
         !(modality %in% c("23","24") & cwt_tx_date < mc$modality_2324_from))
j <- A %>% filter(!is.na(tx_date)) %>%
  left_join(cwt_p %>% select(pseudo_patientid, dtt_date, cwt_tx_date), by = "pseudo_patientid") %>%
  mutate(tx_date_diff = as.numeric(abs(tx_date - cwt_tx_date)),
         wt_dx_to_tx  = as.numeric(tx_date - diagmdy))
funnel <- tibble(
  step = c("left (NCRAS+HES)","linkable CWT record","tx dates agree","non-negative wait","deduplicated"),
  n_patients = c(
    n_distinct(A$pseudo_patientid),
    n_distinct(j$pseudo_patientid[!is.na(j$tx_date_diff)]),
    n_distinct(j$pseudo_patientid[!is.na(j$tx_date_diff) & j$tx_date_diff <= mc$tx_date_diff_max]),
    n_distinct(j$pseudo_patientid[!is.na(j$tx_date_diff) & j$tx_date_diff <= mc$tx_date_diff_max & j$wt_dx_to_tx >= 0]),
    nrow(syn))
)
cat("== Exclusion funnel (from saved tables) ==\n")
print(funnel)