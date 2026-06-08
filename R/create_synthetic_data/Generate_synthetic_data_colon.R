# =============================================================================
# Colon waiting times - synthetic data generator
# -----------------------------------------------------------------------------
# Builds two synthetic tables that mirror the real pipeline inputs:
#   Table A  syn_ncras_hes  - one row per patient (NCRAS registry + the linked
#                             HES resection episode), the LEFT side of the merge
#   Table B  syn_cwt         - raw CWT records (multiple rows per patient,
#                             dd/mm/YYYY date strings), the RIGHT side
#
# It then re-runs the real cancer-waiting-times linkage (the 1b step) on the
# synthetic tables to produce the analysis cohort and the patient-exclusion
# funnel, so the same flow chart can be reproduced off-server.
#
# Driven by colon_profile_for_synthetic.rds when present; otherwise uses the
# built-in defaults so it runs standalone. No real data is used.
#
# Deps: tidyverse, lubridate, haven.
# =============================================================================

library(tidyverse)
library(lubridate)
library(haven)

base_dir <- "Data/synthetic/"
set.seed(2045)

# left-side scale (NCRAS+HES patients entering the CWT linkage); the final
# cohort is smaller after the exclusions below. Override as you like.
N_TOTAL <- 55000L

prof_path <- paste0(base_dir, "colon_profile_for_synthetic.rds")
spec_path <- paste0(base_dir, "colon_pipeline_spec.rds")
profile   <- if (file.exists(prof_path)) readRDS(prof_path) else NULL
spec      <- if (file.exists(spec_path)) readRDS(spec_path) else NULL

colon_icd <- if (!is.null(spec)) spec$colon_icd else
  c("C18","C180","C181","C182","C183","C184","C185","C186","C187","C188","C189")
mc <- if (!is.null(spec)) spec$merge_const else list(
  tx_window_days = c(-90L,180L), tx_date_diff_max = 5L,
  modality_keep = c("01","23","24"), modality_2324_from = as.Date("2020-06-01"))

# -----------------------------------------------------------------------------
# Built-in defaults (used for any profile section that is absent)
# -----------------------------------------------------------------------------
def_marg <- function(levels, props) tibble(level = levels, prop = props / sum(props))
defaults <- list(
  by_year   = def_marg(2015:2022, c(.11,.12,.12,.13,.13,.10,.145,.145)),
  sex       = def_marg(c("Male","Female"), c(.54,.46)),
  age_grp   = def_marg(c("<50","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"),
                       c(.06,.05,.07,.10,.14,.17,.16,.14,.11)),
  ethnicity = def_marg(c("White","Asian","Black","Mixed/Multiple","Other","Unknown"),
                       c(.88,.045,.025,.012,.018,.02)),
  imd       = def_marg(1:5, c(.205,.205,.20,.195,.195)),
  stage     = def_marg(c("1","2","3"), c(.27,.34,.39)),
  route     = def_marg(c("Two Week Wait","GP referral","Screening","Other outpatient","Inpatient elective"),
                       c(.52,.20,.14,.10,.04)),
  screen    = def_marg(c(0,1), c(.86,.14)),
  proc      = def_marg(c("right_hemi","left_hemi","sigmoid","transverse","total_subtotal"),
                       c(.42,.14,.22,.06,.16)),
  cci       = def_marg(c("0","1","2","3+"), c(.66,.18,.10,.06)),
  perf      = def_marg(c("0","1","2","3","4"), c(.55,.22,.13,.07,.03)),
  dco       = def_marg(c(0,1), c(.99,.01)),
  emergency = def_marg(c(0,1), c(.96,.04)),
  modality  = def_marg(c("01","23","24"), c(.90,.06,.04)),
  site_icd  = def_marg(c("C180","C181","C182","C183","C184","C185","C186","C189"),
                       c(.18,.17,.10,.20,.08,.07,.10,.10)),
  mean_dx_to_dtt = 22, sd_dx_to_dtt = 18,
  mean_dtt_to_tx = 24, sd_dtt_to_tx = 20,
  sd_hosp_dtt = 4, sd_hosp_tx = 3,
  n_trust = 110L, n_hosp = 148L,
  change_trust_rate = 0.18,
  cwt_coverage = 0.92, pct_exact = 0.60, pct_within_5 = 0.85,
  pct_dtt_missing = 0.10, pct_mdt = 0.55,
  records_per_patient = def_marg(1:3, c(.82,.14,.04))
)

# pull a marginal from the profile (level/prop), else fall back to default
marg <- function(key_profile, key_default) {
  m <- tryCatch(profile$marginals[[key_profile]], error = function(e) NULL)
  if (is.null(m) || !nrow(m)) defaults[[key_default]] else
    m %>% transmute(level, prop = prop / sum(prop))
}
num_or <- function(x, d) if (is.null(x) || is.na(x)) d else x
interval_ms <- function(key, m_def, s_def) {  # mean+sd of an interval, with fallback
  q <- tryCatch(profile$intervals_overall[[key]], error = function(e) NULL)
  if (is.null(q)) c(m_def, s_def) else c(q$mean, q$sd)
}

# -----------------------------------------------------------------------------
# Samplers
# -----------------------------------------------------------------------------
sample_cat <- function(tbl, n) {
  p <- tbl$prop; p[is.na(p)] <- 0
  sample(tbl$level, n, replace = TRUE, prob = p / sum(p))
}
# gamma draw matched to a target mean and sd (>= 0)
rgamma_ms <- function(n, m, s) {
  m <- max(m, 0.5); s <- max(s, 1)
  shape <- (m / s)^2; scale <- s^2 / m
  pmax(0L, as.integer(round(rgamma(n, shape = shape, scale = scale))))
}
fmt <- function(d) ifelse(is.na(d), NA_character_, format(d, "%d/%m/%Y"))

# =============================================================================
# Trust / hospital pools (NHS ODS style; hospitals nested within trusts)
# =============================================================================
n_trust <- as.integer(num_or(profile$volume$diag_trust$n_distinct, defaults$n_trust))
n_hosp  <- as.integer(num_or(profile$volume$diag_hosp$n_distinct,  defaults$n_hosp))
n_hosp  <- max(n_hosp, n_trust)

chars     <- c(LETTERS, 0:9)
all_2char <- as.vector(outer(chars, chars, paste0))
trust_codes <- paste0("R", sample(all_2char, n_trust))         # e.g. RX9, RJ1

trust_w <- sort(rlnorm(n_trust, 1, 1), decreasing = TRUE); trust_w <- trust_w / sum(trust_w)
hosp_per_trust <- rep(1L, n_trust)
extra <- n_hosp - n_trust
if (extra > 0)
  hosp_per_trust <- hosp_per_trust +
  tabulate(sample(seq_len(n_trust), extra, TRUE, prob = trust_w), nbins = n_trust)

hosp_codes <- unlist(mapply(function(t, k) paste0(t, sample(all_2char, k)),
                            trust_codes, hosp_per_trust, SIMPLIFY = FALSE))
hosp_trust <- rep(trust_codes, times = hosp_per_trust)
hosp_w     <- rlnorm(length(hosp_codes), 0, 0.8)
hosp_w     <- ave(hosp_w, hosp_trust, FUN = function(x) x / sum(x))

draw_trust <- function(n) sample(trust_codes, n, TRUE, trust_w)
draw_hosp  <- function(trusts) {
  out <- character(length(trusts))
  for (t in unique(trusts)) {
    idx <- which(trusts == t)
    h <- hosp_codes[hosp_trust == t]; w <- hosp_w[hosp_trust == t]
    out[idx] <- sample(h, length(idx), TRUE, w)
  }
  out
}
# per-hospital random intercepts on the two wait components (the ICC signal)
re_dtt <- setNames(rnorm(length(hosp_codes), 0, num_or(profile$between_hosp_sd$wt_dx_to_dtt, defaults$sd_hosp_dtt)), hosp_codes)
re_tx  <- setNames(rnorm(length(hosp_codes), 0, num_or(profile$between_hosp_sd$wt_dtt_to_tx, defaults$sd_hosp_tx)), hosp_codes)

# =============================================================================
# Table A: NCRAS + HES cohort (left side of the CWT merge)
# =============================================================================
n <- N_TOTAL
A <- tibble(
  pseudo_patientid = sprintf("S%07d", seq_len(n)),
  pseudo_tumourid  = sprintf("T%07d", seq_len(n))
)

yr_tbl <- if (!is.null(profile$by_year)) profile$by_year %>% transmute(level, prop = prop / sum(prop)) else defaults$by_year
A$ydiag <- as.integer(as.character(sample_cat(yr_tbl, n)))
A$diagmdy <- as.Date(paste0(A$ydiag, "-01-01")) + as.integer(runif(n, 0, 364))

# age: band then a value inside it
age_grp <- sample_cat(marg("age_grp","age_grp"), n)
lo <- c("<50"=40,"50-54"=50,"55-59"=55,"60-64"=60,"65-69"=65,"70-74"=70,"75-79"=75,"80-84"=80,"85+"=85)
hi <- c("<50"=49,"50-54"=54,"55-59"=59,"60-64"=64,"65-69"=69,"70-74"=74,"75-79"=79,"80-84"=84,"85+"=95)
A$agediag <- as.integer(round(runif(n, lo[age_grp], hi[age_grp])))

A$sex                              <- sample_cat(marg("sex","sex"), n)
A$ethnicity_group_broad            <- sample_cat(marg("ethnicity","ethnicity"), n)

# imd: profile may store "1 - most deprived" / "5 - least deprived" as labels;
# extract just the leading digit so as.integer() is reliable
A$NHSE_reversed_imd_quintile_lsoas <- as.integer(substr(
  as.character(sample_cat(marg("imd_quintile","imd"), n)), 1, 1
))

A$stage                            <- as.character(sample_cat(marg("stage","stage"), n))

# route: profile may store "TWW" rather than "Two Week Wait"
A$route_combined <- recode(
  as.character(sample_cat(marg("route_combined","route"), n)),
  "TWW" = "Two Week Wait"
)

# screendetected: profile stores "Y"/"N"/"" — map to 1/0 explicitly
scr_raw          <- as.character(sample_cat(marg("screendetected","screen"), n))
A$screendetected <- case_when(scr_raw == "Y" ~ 1L, scr_raw == "N" ~ 0L, TRUE ~ 0L)

A$colon_proc_type_primary          <- sample_cat(marg("proc_type","proc"), n)
A$cci_group                        <- as.character(sample_cat(marg("cci_group","cci"), n))
A$perf_status                      <- as.integer(sample_cat(marg("perf_status","perf"), n))

# dco: profile may store "Y"/"N" — map to 1/0 explicitly
dco_raw  <- as.character(sample_cat(marg("dco","dco"), n))
A$dco    <- case_when(dco_raw == "Y" ~ 1L, dco_raw == "N" ~ 0L, TRUE ~ 0L)

A$emergency                        <- as.integer(sample_cat(marg("emergency","emergency"), n))

# yearmonth_diag: string "YYYYMM", used by the Stata analysis_year derivation
A$yearmonth_diag <- format(A$diagmdy, "%Y%m")

# metastasis exclusion stubs - cohort is stage 1-3 / M0 by construction, so
# all are 0, but the Stata "drop if any_mets==1 / tnm_m==1" lines need them
A$any_mets   <- 0L
A$tnm_m      <- 0L
A$pretreat_m <- 0L

# rcs_ch_score: numeric Charlson score (the rapid analysis uses i.rcs_ch_score)
A$rcs_ch_score <- case_when(
  A$cci_group == "0"  ~ 0L,
  A$cci_group == "1"  ~ 1L,
  A$cci_group == "2"  ~ 2L,
  A$cci_group == "3+" ~ 3L,
  TRUE                ~ NA_integer_
)

# screen-detected mostly comes via the Screening route
A$screendetected <- if_else(A$route_combined == "Screening", 1L, A$screendetected)

# tumour site (drives site_match later); substage label
A$sitestr    <- sample_cat(defaults$site_icd, n)
A$stage_best <- paste0(A$stage, sample(c("","A","B"), n, TRUE, c(.6,.25,.15)))

# trust / hospital, nested; minority diagnosed and treated at different trusts
A$diag_trust  <- draw_trust(n)
A$first_trust <- A$diag_trust
A$diag_hosp   <- draw_hosp(A$diag_trust)
moved         <- runif(n) < num_or(profile$change_trust_rate, defaults$change_trust_rate)
treat_trust   <- A$diag_trust
treat_trust[moved] <- draw_trust(sum(moved))
A$PROCODE3    <- treat_trust
A$SITETRET    <- if_else(moved, draw_hosp(treat_trust), A$diag_hosp)
A$change_trust <- substr(A$diag_trust, 1, 3) != substr(A$PROCODE3, 1, 3)
A$tx_trust     <- A$PROCODE3   # treatment trust, 3-char (from HES PROCODE3)
A$tx_hosp      <- A$SITETRET   # treatment hospital/site, 5-char (from HES SITETRET)

# survival / death flag
ms_surv <- interval_ms("surv_from_dx_days", 900, 600)
A$finmdy <- A$diagmdy + pmax(rgamma_ms(n, ms_surv[1], ms_surv[2]), 1L)
A$dead   <- as.integer(runif(n) < 0.28)
A$finmdy[A$dead == 0] <- as.Date(NA)

# latent waiting-time components (centred covariate effects keep the means on
# target while giving the analysis a signal to recover)
ms_d <- interval_ms("wt_dx_to_dtt", defaults$mean_dx_to_dtt, defaults$sd_dx_to_dtt)
ms_t <- interval_ms("wt_dtt_to_tx", defaults$mean_dtt_to_tx, defaults$sd_dtt_to_tx)
eff_route <- c("Two Week Wait"=0,"GP referral"=8,"Screening"=2,"Other outpatient"=4,"Inpatient elective"=3)
eff_cci   <- c("0"=0,"1"=1,"2"=2,"3+"=3)
eff_stage <- c("1"=0,"2"=-1,"3"=-1.5)
cov_eff <- coalesce(eff_route[A$route_combined], 0) +
  coalesce(eff_cci[A$cci_group], 0) +
  eff_stage[A$stage] +
  0.03 * pmax(0, A$agediag - 70) + (A$ydiag - 2019) * 0.8
cov_eff <- cov_eff - mean(cov_eff, na.rm = TRUE)

dx_to_dtt <- pmax(0L, rgamma_ms(n, ms_d[1], ms_d[2]) + round(re_dtt[A$diag_hosp] + 0.45 * cov_eff))
dtt_to_tx <- pmax(0L, rgamma_ms(n, ms_t[1], ms_t[2]) + round(re_tx[A$diag_hosp]  + 0.55 * cov_eff))
total     <- dx_to_dtt + dtt_to_tx
over      <- which(total > 180)            # keep total wait inside the cohort window
while (length(over)) {
  dx_to_dtt[over] <- pmax(0L, rgamma_ms(length(over), ms_d[1], ms_d[2]) + round(re_dtt[A$diag_hosp[over]]))
  dtt_to_tx[over] <- pmax(0L, rgamma_ms(length(over), ms_t[1], ms_t[2]) + round(re_tx[A$diag_hosp[over]]))
  total[over] <- dx_to_dtt[over] + dtt_to_tx[over]
  over <- over[total[over] > 180]
}

A$tx_date           <- A$diagmdy + total
A$days_diag_to_surg <- as.integer(A$tx_date - A$diagmdy)   # in [0,180]

# carry the latent dtt offset internally to anchor the CWT records (not exported)
lat_dx_to_dtt <- dx_to_dtt

# =============================================================================
# Table B: raw CWT records (right side), anchored to Table A
# =============================================================================
cov     <- num_or(profile$cwt_coverage$pct_any_cwt, defaults$cwt_coverage)
p_exact <- num_or(profile$cwt_agreement$pct_exact,  defaults$pct_exact)
p_w5    <- num_or(profile$cwt_agreement$pct_within_5, defaults$pct_within_5)
p_dttNA <- num_or(profile$pct_dtt_missing_in_cohort, defaults$pct_dtt_missing)
p_mdt   <- num_or(profile$cwt_completeness$pct_mdt, defaults$pct_mdt)
mod_tbl <- if (!is.null(profile$cwt_modality)) profile$cwt_modality %>% transmute(level, prop = prop / sum(prop)) else defaults$modality
recs    <- if (!is.null(profile$cwt_records_per_patient))
  profile$cwt_records_per_patient %>% transmute(level = records, prop) else defaults$records_per_patient

has_cwt <- runif(n) < cov
idx     <- which(has_cwt)

# agreement offset between the CWT treatment date and the HES surgery date
u <- runif(length(idx))
delta <- integer(length(idx))
mid <- u >= p_exact & u < p_w5
far <- u >= p_w5
delta[mid] <- sample(c(-5:-1, 1:5), sum(mid), TRUE)
delta[far] <- sample(c(-30:-6, 6:30), sum(far), TRUE)   # excluded by tx_date_diff > 5

dtt_anchor   <- A$diagmdy[idx] + lat_dx_to_dtt[idx]
treat_anchor <- A$tx_date[idx] + delta
dtt_missing  <- runif(length(idx)) < p_dttNA
dtt_anchor[dtt_missing] <- as.Date(NA)

mdt_have <- runif(length(idx)) < p_mdt
mdt_anchor <- as.Date(rep(NA, length(idx)), origin = "1970-01-01")
mdt_anchor[mdt_have] <- dtt_anchor[mdt_have] - sample(0:21, sum(mdt_have), TRUE)

crtp  <- dtt_anchor - sample(20:60, length(idx), TRUE)
fseen <- crtp + sample(0:14, length(idx), TRUE)
site  <- ifelse(runif(length(idx)) < 0.95,           # site usually matches the registry
                A$sitestr[idx],
                sample(defaults$site_icd$level, length(idx), TRUE))
modal <- sample_cat(mod_tbl, length(idx))

# org_treat_start: the CWT treating organisation code. Mostly matches the HES
# SITETRET (same event, different source), but ~5% differ due to data quality —
# so the tx_trust vs tx_trust_cwt check in the validation is non-trivial.
cwt_org <- A$SITETRET[idx]
noise   <- runif(length(idx)) < 0.05
cwt_org[noise] <- sample(hosp_codes, sum(noise), replace = TRUE)

anchor <- tibble(
  pseudo_patientid   = A$pseudo_patientid[idx],
  site_icd10         = site,
  modality           = as.character(modal),
  org_treat_start    = cwt_org,
  crtp_date          = fmt(crtp),
  date_first_seen    = fmt(fseen),
  mdt_date           = fmt(mdt_anchor),
  treat_period_start = fmt(dtt_anchor),
  treat_start        = fmt(treat_anchor)
)

# a minority of patients carry extra (non-anchor) records placed well away from
# the surgery date, so they are dropped by the tx_date_diff filter and the
# dedup step keeps the anchor
k <- as.integer(sample_cat(recs, length(idx)))
extra_idx <- which(k > 1)
extra <- map_dfr(extra_idx, function(j) {
  m <- k[j] - 1L
  base_dtt <- A$diagmdy[idx[j]] + lat_dx_to_dtt[idx[j]] + cumsum(sample(30:120, m, TRUE))
  tibble(
    pseudo_patientid   = A$pseudo_patientid[idx[j]],
    site_icd10         = anchor$site_icd10[j],
    modality           = as.character(sample_cat(mod_tbl, m)),
    org_treat_start    = NA_character_,
    crtp_date          = fmt(base_dtt - sample(20:60, m, TRUE)),
    date_first_seen    = NA_character_,
    mdt_date           = NA_character_,
    treat_period_start = fmt(base_dtt),
    treat_start        = fmt(base_dtt + sample(10:40, m, TRUE))
  )
})
syn_cwt <- bind_rows(anchor, extra) %>% arrange(pseudo_patientid)

# =============================================================================
# Re-run the real CWT linkage (1b) and record the exclusion funnel
# =============================================================================
run_merge <- function(A, cwt) {
  cwt_p <- cwt %>%
    mutate(dtt_date    = as.Date(treat_period_start, "%d/%m/%Y"),
           cwt_tx_date = as.Date(treat_start,        "%d/%m/%Y"),
           mdt_date    = as.Date(mdt_date,           "%d/%m/%Y")) %>%
    filter(modality %in% mc$modality_keep,
           !(modality %in% c("23","24") & cwt_tx_date < mc$modality_2324_from))
  
  step0 <- A %>% filter(!is.na(tx_date))
  step1 <- step0 %>%
    left_join(cwt_p, by = "pseudo_patientid") %>%
    mutate(tx_date_diff = as.numeric(abs(tx_date - cwt_tx_date)))
  step2 <- step1 %>% filter(!is.na(tx_date_diff))                 # a linkable record
  step3 <- step2 %>% filter(tx_date_diff <= mc$tx_date_diff_max)  # date agreement
  step4 <- step3 %>%
    mutate(
      dx_le_mdt_ok  = is.na(mdt_date) | diagmdy <= mdt_date,
      dx_le_dtt_ok  = is.na(dtt_date) | diagmdy <= dtt_date,
      dx_le_tx_ok   = diagmdy <= tx_date,
      mdt_le_dtt_ok = is.na(mdt_date) | is.na(dtt_date) | mdt_date <= dtt_date,
      dtt_le_tx_ok  = is.na(dtt_date) | dtt_date <= tx_date,
      mdt_le_tx_ok  = is.na(mdt_date) | mdt_date <= tx_date,
      seq_ok        = dx_le_dtt_ok & dx_le_mdt_ok & dx_le_tx_ok &
        mdt_le_dtt_ok & dtt_le_tx_ok & mdt_le_tx_ok,
      wt_dx_to_dtt  = as.numeric(dtt_date - diagmdy),
      wt_dx_to_tx   = as.numeric(tx_date  - diagmdy),
      wt_dtt_to_tx  = as.numeric(tx_date  - dtt_date),
      site_match    = as.integer(str_sub(site_icd10, 1, 3) == str_sub(sitestr, 1, 3)),
      tx_trust_cwt  = substr(org_treat_start, 1, 3),   # treatment trust from CWT source
      # signed date difference (mirrors Stata's diff_cwt_cr_treat = date_of_surgery - date_treat_start)
      diff_cwt_cr_treat = as.integer(tx_date - cwt_tx_date),
      diff_cwt_cr_treat_cat = case_when(
        is.na(cwt_tx_date)                                             ~ 5L,
        diff_cwt_cr_treat == 0                                         ~ 1L,
        diff_cwt_cr_treat >  0 & diff_cwt_cr_treat <  5               ~ 2L,
        diff_cwt_cr_treat < -0 & diff_cwt_cr_treat > -5               ~ 2L,
        diff_cwt_cr_treat >  5 & diff_cwt_cr_treat <  30              ~ 3L,
        diff_cwt_cr_treat < -5 & diff_cwt_cr_treat > -30              ~ 3L,
        diff_cwt_cr_treat >  30 | diff_cwt_cr_treat < -30             ~ 4L,
        TRUE                                                           ~ 5L
      )
    ) %>%
    filter(wt_dx_to_tx >= 0)
  final <- step4 %>%
    arrange(pseudo_patientid, desc(site_match), tx_date_diff, dtt_date, cwt_tx_date) %>%
    distinct(pseudo_patientid, .keep_all = TRUE) %>%
    mutate(missing_CWT = 0L)   # all retained patients have a CWT record
  
  funnel <- tibble(
    step = c("NCRAS+HES patients (left)",
             "with a linkable CWT record (modality kept)",
             "treatment dates agree (tx_date_diff <= 5)",
             "non-negative total wait (wt_dx_to_tx >= 0)",
             "deduplicated to one record per patient"),
    n_patients = c(n_distinct(step0$pseudo_patientid),
                   n_distinct(step2$pseudo_patientid),
                   n_distinct(step3$pseudo_patientid),
                   n_distinct(step4$pseudo_patientid),
                   n_distinct(final$pseudo_patientid)),
    n_rows = c(nrow(step0), nrow(step2), nrow(step3), nrow(step4), nrow(final))
  )
  list(cohort = final, funnel = funnel)
}

merged <- run_merge(A, syn_cwt)
syn_cohort <- merged$cohort

cat("\nCWT linkage funnel (synthetic):\n")
print(merged$funnel)

# =============================================================================
# Select the analysis columns and save
# =============================================================================
keep <- c(
  "pseudo_patientid","pseudo_tumourid","diagmdy","ydiag","yearmonth_diag",
  "sitestr","stage","stage_best",
  "agediag","sex","ethnicity_group_broad","NHSE_reversed_imd_quintile_lsoas",
  "route_combined","screendetected","cci_group","rcs_ch_score","perf_status",
  "colon_proc_type_primary","any_mets","tnm_m","pretreat_m",
  "diag_hosp","diag_trust","first_trust","tx_date","days_diag_to_surg",
  "tx_hosp","tx_trust","SITETRET","PROCODE3",
  "emergency","dco","dead","finmdy","change_trust",
  "site_icd10","modality","dtt_date","cwt_tx_date","mdt_date",
  "tx_date_diff","diff_cwt_cr_treat","diff_cwt_cr_treat_cat","missing_CWT","site_match",
  "tx_trust_cwt",
  "dx_le_mdt_ok","dx_le_dtt_ok","dx_le_tx_ok","mdt_le_dtt_ok","dtt_le_tx_ok","mdt_le_tx_ok",
  "seq_ok","wt_dx_to_dtt","wt_dtt_to_tx","wt_dx_to_tx"
)
syn_cohort <- syn_cohort %>% select(any_of(keep))

A_out <- A   # full left-side table, lets you re-run the merge yourself

saveRDS(A_out,       paste0(base_dir, "colon_ncras_hes_synthetic.rds"))
saveRDS(syn_cwt,     paste0(base_dir, "colon_cwt_records_synthetic.rds"))
saveRDS(syn_cohort,  paste0(base_dir, "colon_cohort_synthetic.rds"))

to_stata <- function(df) df %>% mutate(across(where(is.factor), as.character),
                                       across(where(is.logical), as.integer))
write_dta(to_stata(A_out),      paste0(base_dir, "colon_ncras_hes_synthetic.dta"))
write_dta(to_stata(syn_cwt),    paste0(base_dir, "colon_cwt_records_synthetic.dta"))
write_dta(to_stata(syn_cohort), paste0(base_dir, "colon_cohort_synthetic.dta"))

cat("\nSaved Table A (", nrow(A_out), "rows), Table B CWT (", nrow(syn_cwt),
    "rows), analysis cohort (", nrow(syn_cohort), "rows).\n")
cat("Means dx_to_dtt / dtt_to_tx / dx_to_tx:",
    round(mean(syn_cohort$wt_dx_to_dtt, na.rm = TRUE), 1),
    round(mean(syn_cohort$wt_dtt_to_tx, na.rm = TRUE), 1),
    round(mean(syn_cohort$wt_dx_to_tx), 1), "\n")