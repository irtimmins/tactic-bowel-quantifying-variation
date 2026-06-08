# =============================================================================
# 01_extract_icon_metadata.R
# -----------------------------------------------------------------------------
# RUN THIS ON THE SECURE ICON SERVER, against the real colon cohort.
#
# It does TWO things and writes only AGGREGATE, NON-DISCLOSIVE output:
#   1. Validates that the real cohort has the expected 172-column structure and
#      writes a structural data dictionary (classes + rounded counts only - NO
#      values) to `column_spec.csv`.
#   2. Extracts summary parameters (category proportions, waiting-time moments,
#      between-hospital variance, year/covariate effects) into `icon_metadata.rds`
#      and a human-readable `icon_metadata_summary.txt`.
#
# Small cells are suppressed (counts < MIN_CELL are merged/dropped) and all
# counts are rounded, so the outputs are safe to take through the standard
# disclosure check and off the server. NO row-level data is ever written.
#
# The resulting icon_metadata.rds is then consumed by 02_build_synthetic_cohort.R
# Dependencies: base R only.
# =============================================================================

# ----------------------------- USER SETTINGS ---------------------------------
COHORT_PATH <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_2015_2022.rds"
CCI_PATH    <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_cci_2015_2022.rds" # optional
OUT_DIR     <- "."
MIN_CELL    <- 10L     # suppress / merge categories with fewer than this many obs
ROUND_TO    <- 10L     # round all disclosed counts to nearest this
# -----------------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
d <- readRDS(COHORT_PATH)
cci <- if (file.exists(CCI_PATH)) readRDS(CCI_PATH) else NULL

# --------------------- 1. STRUCTURE VALIDATION + SPEC ------------------------
expected_172 <- c(
 "pseudo_patientid","pseudo_tumourid","diagmdy","ydiag","cancer","sitestr","typestr",
 "basisofdiagnosis","grade","stage_best","stage_best_system","t_best","n_best","m_best",
 "t_path","n_path","m_path","sex","agediag","birthmdy","ethnicity_group_broad","lsoa11_code",
 "NHSE_reversed_imd_quintile_lsoas","canalliance_2024_code","canalliance_2024_name",
 "diag_trust","diag_trust_name","first_trust","first_trust_name","first_hosp_date",
 "diag_hosp","diag_hosp_name","route_bjc","final_route","route_code","tww_to_treat",
 "sg_flag","rt_flag","ct_flag","screendetected","dead","finmdy","dco","er_status",
 "pr_status","her2_status","laterality.x","dukes","nodesexcised","nodesinvolved",
 "final_route_chr","route_bjc_chr","route_combined","stage","STUDY_ID","ADMIDATE",
 "ADMIMETH","PROCODE3","SITETRET","EPISTART","EPIORDER","EPITYPE",
 sprintf("OPDATE_%02d", 1:24),
 "op_position","colon_opcs_primary","opcs3","colon_opcs_primary_position",
 "colon_proc_type_primary","primary_flag","all_colon_opcs","all_colon_proc_types",
 "all_colon_opcs_fields","n_colon_codes_in_episode","emergency","days_diag_to_surg",
 "tx_date","org_ppi","ref_source","priority_type","dec_to_ref_date","month_dec_to_ref_date",
 "year_dec_to_ref_date","crtp_date","month_crtp_date","year_crtp_date","ref_type",
 "cons_upgrade_date","month_cons_upgrade_date","year_cons_upgrade_date","org_cons_upgrade",
 "delay_cons_treat_reason","date_first_seen","month_date_first_seen","year_date_first_seen",
 "org_first_seen","wta_first_seen","wta_first_seen_reason","delay_ref_fs_reason",
 "patient_status","laterality.y","mets_site","treat_period_start","month_treat_period_start",
 "year_treat_period_start","org_dec_to_treat","fdp_end_reason","fdp_diag_site","fdp_end_date",
 "month_fdp_end_date","year_fdp_end_date","delay_fdp_reason","fdp_exclusion_reason",
 "fdp_outcome_prof_type","fdp_outcome_method","org_fdp_end","treat_start","month_treat_start",
 "year_treat_start","org_treat_start","cte_type","modality","clin_trial","care_setting",
 "delay_dtt_treat_reason","wta_treat","wta_treat_reason","delay_ref_treat_reason",
 "radio_intent","radio_priority","mdt_ind","mdt_date","month_mdt_date","year_mdt_date",
 "practice_code","site_icd10","dtt_date","cwt_tx_date","tx_date_diff","dx_le_mdt_ok",
 "dx_le_dtt_ok","dx_le_tx_ok","mdt_le_dtt_ok","dtt_le_tx_ok","mdt_le_tx_ok","seq_ok",
 "wt_dx_to_dtt","wt_dx_to_tx","wt_dtt_to_tx","site_match")

cat("== STRUCTURE CHECK ==\n")
cat("cols in cohort:", ncol(d), " | expected:", length(expected_172), "\n")
missing_cols <- setdiff(expected_172, names(d))
extra_cols   <- setdiff(names(d), expected_172)
if (length(missing_cols)) cat("MISSING vs expected:", paste(missing_cols, collapse=", "), "\n")
if (length(extra_cols))   cat("EXTRA   vs expected:", paste(extra_cols,   collapse=", "), "\n")

rnd <- function(x) round(x / ROUND_TO) * ROUND_TO
spec <- data.frame(
  column   = names(d),
  class    = vapply(d, function(x) class(x)[1], character(1)),
  n_missing_rounded = vapply(d, function(x) rnd(sum(is.na(x))), numeric(1)),
  n_unique_rounded  = vapply(d, function(x) rnd(length(unique(x))), numeric(1)),
  in_expected_172   = names(d) %in% expected_172,
  row.names = NULL, stringsAsFactors = FALSE
)
write.csv(spec, file.path(OUT_DIR, "column_spec.csv"), row.names = FALSE)
cat("Wrote column_spec.csv (", nrow(spec), "rows)\n\n")

# ----------------------- 2. HELPER: SAFE PROPORTIONS -------------------------
# category proportions with small-cell suppression (cells < MIN_CELL dropped,
# remainder renormalised, proportions rounded to 3 dp)
safe_props <- function(x, keep_levels = NULL) {
  x <- as.character(x); x <- x[!is.na(x) & x != ""]
  tb <- table(x)
  tb <- tb[tb >= MIN_CELL]
  if (!length(tb)) return(NULL)
  if (!is.null(keep_levels)) tb <- tb[names(tb) %in% keep_levels]
  p <- as.numeric(tb) / sum(tb)
  setNames(round(p, 3), names(tb))
}
# mean difference of a numeric outcome by group vs a reference level (rounded)
eff_vs_ref <- function(num, grp, ref) {
  grp <- as.character(grp); ok <- !is.na(num) & !is.na(grp)
  num <- num[ok]; grp <- grp[ok]
  m <- tapply(num, grp, function(z) if (length(z) >= MIN_CELL) mean(z) else NA_real_)
  if (!(ref %in% names(m)) || is.na(m[[ref]])) return(NULL)
  round(m - m[[ref]], 2)
}
# one-way between-group SD via ANOVA method of moments (random-intercept proxy)
between_sd <- function(num, grp) {
  ok <- !is.na(num) & !is.na(grp); num <- num[ok]; grp <- as.character(grp)[ok]
  if (length(unique(grp)) < 2) return(NA_real_)
  fit <- summary(aov(num ~ grp))[[1]]
  msb <- fit["grp", "Mean Sq"]; msw <- fit["Residuals", "Mean Sq"]
  ni  <- table(grp); n0 <- (sum(ni) - sum(ni^2)/sum(ni)) / (length(ni) - 1)
  sqrt(max(0, (msb - msw) / n0))
}

# ----------------------------- 3. PARAMETERS ---------------------------------
yrs   <- as.integer(d$ydiag)
yr_tb <- table(yrs[!is.na(yrs)])
mom_shape <- function(x) { x <- x[!is.na(x) & x >= 0]; m <- mean(x); v <- var(x); max(1.1, m^2 / v) }

params <- list(
  year_levels  = as.integer(names(yr_tb)),
  year_weights = round(as.numeric(yr_tb)/sum(yr_tb), 4),

  mean_dx_to_dtt = round(mean(d$wt_dx_to_dtt[d$wt_dx_to_dtt >= 0], na.rm = TRUE), 2),
  mean_dtt_to_tx = round(mean(d$wt_dtt_to_tx[d$wt_dtt_to_tx >= 0], na.rm = TRUE), 2),
  gamma_shape_dtt = round(mom_shape(d$wt_dx_to_dtt), 2),
  gamma_shape_tx  = round(mom_shape(d$wt_dtt_to_tx), 2),
  sd_hosp_dtt = round(between_sd(d$wt_dx_to_dtt, d$diag_hosp), 2),
  sd_hosp_tx  = round(between_sd(d$wt_dtt_to_tx, d$diag_hosp), 2),

  year_effect = eff_vs_ref(d$wt_dx_to_tx, yrs, ref = "2019"),

  p_sex        = safe_props(d$sex),
  p_ethnicity  = safe_props(d$ethnicity_group_broad),
  p_imd        = safe_props(d$NHSE_reversed_imd_quintile_lsoas),
  p_stage      = safe_props(d$stage),
  p_route      = safe_props(d$route_combined),
  p_grade      = safe_props(d$grade),
  p_proc       = safe_props(d$colon_proc_type_primary),
  p_modality   = safe_props(d$modality),

  eff_route    = eff_vs_ref(d$wt_dx_to_tx, d$route_combined, ref = "Two Week Wait"),
  eff_stage    = eff_vs_ref(d$wt_dx_to_tx, d$stage,          ref = "1"),
  age_slope    = round(unname(coef(lm(wt_dx_to_tx ~ I(pmax(0, agediag - 70)), data = d))[2]), 3),

  p_emergency  = round(mean(d$emergency, na.rm = TRUE), 3),
  p_dead       = round(mean(as.numeric(d$dead), na.rm = TRUE), 3),
  p_screen_given_route_screening =
    round(mean(as.numeric(d$screendetected)[d$route_combined == "Screening"], na.rm = TRUE), 3),
  mean_age = round(mean(d$agediag, na.rm = TRUE), 1),
  sd_age   = round(sd(d$agediag,   na.rm = TRUE), 1)
)

# CCI distribution / effect from the optional cci file
if (!is.null(cci) && "cci_group" %in% names(cci)) {
  params$p_cci <- safe_props(cci$cci_group)
  if ("wt_dx_to_tx" %in% names(cci))
    params$eff_cci <- eff_vs_ref(cci$wt_dx_to_tx, cci$cci_group, ref = "0")
}

meta <- list(
  n_patients  = rnd(nrow(d)),
  n_hospitals = rnd(length(unique(d$diag_hosp[grepl("^R[A-Z0-9]{4}$", d$diag_hosp)]))),
  n_trusts    = rnd(length(unique(d$diag_trust))),
  generated   = as.character(Sys.Date()),
  min_cell    = MIN_CELL,
  params      = params
)

saveRDS(meta, file.path(OUT_DIR, "icon_metadata.rds"))

# ----------------------- 4. HUMAN-READABLE SUMMARY ---------------------------
sink(file.path(OUT_DIR, "icon_metadata_summary.txt"))
cat("ICON colon cohort - aggregate metadata for synthetic data generation\n")
cat("Generated:", meta$generated, "| min cell:", MIN_CELL, "| counts rounded to", ROUND_TO, "\n")
cat(strrep("-", 70), "\n")
cat("n patients ~", meta$n_patients, "| n hospitals ~", meta$n_hospitals,
    "| n trusts ~", meta$n_trusts, "\n\n")
str(meta$params)
sink()

cat("\nWrote icon_metadata.rds and icon_metadata_summary.txt\n")
cat("Review column_spec.csv + icon_metadata_summary.txt for disclosure, then\n")
cat("copy icon_metadata.rds out and point 02_build_synthetic_cohort.R at it.\n")
