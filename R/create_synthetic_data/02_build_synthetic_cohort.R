# =============================================================================
# 02_build_synthetic_cohort.R
# -----------------------------------------------------------------------------
# Build a SYNTHETIC colon-cancer waiting-times cohort with the EXACT 172-column
# structure of the real ICON `colon_cohort_2015_2022.rds` (see scripts 1a + 1b).
#
# NO real patient data is used or required to run this script. It can be run:
#   (a) STANDALONE, using the built-in default parameters (derived from the
#       published Stata logs / provider summaries); or
#   (b) DATA-INFORMED, by pointing `META_PATH` at the aggregate
#       `icon_metadata.rds` produced by 01_extract_icon_metadata.R on the
#       secure server (no row-level data leaves the server - only summaries).
#
# Output: synthetic_colon_cohort.dta  (+ .rds)
#
# Dependencies: base R only. `haven` is used for .dta if installed; otherwise
# falls back to `foreign::write.dta`.
# =============================================================================

# ----------------------------- USER SETTINGS ---------------------------------
OUT_DIR    <- "."                       # where to write outputs
META_PATH  <- NULL                      # e.g. "icon_metadata.rds"; NULL = defaults
N_PATIENTS <- 50000L                    # cohort size (real cohort ~55-67k)
N_HOSPITALS<- 148L                      # distinct diagnosing hospitals (~135-150)
SEED       <- 2045L
# -----------------------------------------------------------------------------

set.seed(SEED)

# ------------------------- 0. DEFAULT PARAMETERS -----------------------------
# These are sensible, NON-disclosive defaults taken from the analysis logs.
# If META_PATH is supplied they are overwritten by the real aggregate metadata.
defaults <- list(
  year_levels  = 2015:2022,
  year_weights = c(0.11,0.12,0.12,0.13,0.13,0.10,0.145,0.145),  # dip in 2020 (COVID)

  # Waiting-time component means (days). wt_dx_to_tx = dx_to_dtt + dtt_to_tx.
  mean_dx_to_dtt = 22,   # diagnosis -> decision to treat
  mean_dtt_to_tx = 24,   # decision to treat -> surgery
  gamma_shape_dtt = 2.2, # shape of Gamma for each component (controls skew)
  gamma_shape_tx  = 2.6,
  sd_hosp_dtt = 4.0,     # between-hospital SD (days) on the dx->dtt component
  sd_hosp_tx  = 3.0,     # between-hospital SD (days) on the dtt->tx component

  # Calendar-time effect on total wait (added days vs 2019 reference): COVID era longer
  year_effect = c(`2015`=-1.3,`2016`=-1.2,`2017`=-1.1,`2018`=-0.9,`2019`=0,
                  `2020`=3.4,`2021`=5.4,`2022`=5.8),

  # Categorical distributions (probabilities)
  p_sex        = c(Male=0.54, Female=0.46),
  p_ethnicity  = c(White=0.88, Asian=0.045, Black=0.025, `Mixed/Multiple`=0.012,
                   Other=0.018, Unknown=0.02),
  p_imd        = c(`1`=0.205,`2`=0.205,`3`=0.20,`4`=0.195,`5`=0.195), # 1=most deprived
  p_stage      = c(`1`=0.27,`2`=0.34,`3`=0.39),
  p_route      = c(`Two Week Wait`=0.52,`GP referral`=0.20,`Screening`=0.14,
                   `Other outpatient`=0.10,`Inpatient elective`=0.04),
  p_cci        = c(`0`=0.66,`1`=0.18,`2`=0.10,`3+`=0.06),
  p_grade      = c(G1=0.10,G2=0.62,G3=0.20,GX=0.08),
  p_proc       = c(right_hemi=0.42,left_hemi=0.14,sigmoid=0.22,transverse=0.06,
                   total_subtotal=0.16),
  p_modality   = c(`01`=0.90,`23`=0.06,`24`=0.04),

  # Covariate effects on TOTAL wait (days) - modest, for sign-recovery testing
  eff_route    = c(`Two Week Wait`=0,`GP referral`=8,`Screening`=2,
                   `Other outpatient`=4,`Inpatient elective`=3),
  eff_cci      = c(`0`=0,`1`=1,`2`=2,`3+`=3),
  eff_stage    = c(`1`=0,`2`=-1.0,`3`=-1.5),
  age_slope    = 0.03,    # extra days per year of age above 70
  p_emergency  = 0.04,    # admimeth emergency (cohort is mostly elective)
  p_dead       = 0.28,    # any death during follow-up
  p_screen_given_route_screening = 0.97,
  mean_age = 71, sd_age = 11
)

# Merge in real metadata if provided -----------------------------------------
P <- defaults
if (!is.null(META_PATH) && file.exists(META_PATH)) {
  meta <- readRDS(META_PATH)
  for (nm in intersect(names(meta$params), names(P))) P[[nm]] <- meta$params[[nm]]
  if (!is.null(meta$n_patients))  N_PATIENTS  <- meta$n_patients
  if (!is.null(meta$n_hospitals)) N_HOSPITALS <- meta$n_hospitals
  message("Loaded aggregate metadata from ", META_PATH)
} else {
  message("No metadata supplied - using built-in default parameters.")
}

n <- N_PATIENTS
psamp <- function(probs, size) sample(names(probs), size, TRUE, prob = unname(probs))

# --------------------- 1. HOSPITAL / TRUST STRUCTURE -------------------------
# 3-char trust codes -> 5-char hospital site codes (^R[A-Z0-9]{4}$).
mk_code3 <- function(k) {
  L <- LETTERS; D <- c(LETTERS, 0:9)
  paste0("R", L[((seq_len(k)-1) %% 26) + 1], D[((seq_len(k)*7-1) %% 36) + 1])
}
n_trust <- max(60L, round(N_HOSPITALS * 0.75))
trust_codes <- unique(mk_code3(n_trust))
n_trust <- length(trust_codes)
# assign each hospital to a trust (some trusts have >1 site)
hosp_trust <- sample(trust_codes, N_HOSPITALS, replace = TRUE)
hosp_suffix <- sprintf("%02d", ave(seq_len(N_HOSPITALS), hosp_trust, FUN = seq_along))
hosp_codes  <- paste0(hosp_trust, hosp_suffix)         # e.g. "RBA01"
trust_name  <- setNames(paste0("NHS Trust ", match(trust_codes, trust_codes)), trust_codes)
hosp_name   <- setNames(paste0("Hospital Site ", seq_along(hosp_codes)), hosp_codes)

# hospital random effects (between-hospital variation in the two wait components)
hosp_re_dtt <- setNames(rnorm(N_HOSPITALS, 0, P$sd_hosp_dtt), hosp_codes)
hosp_re_tx  <- setNames(rnorm(N_HOSPITALS, 0, P$sd_hosp_tx),  hosp_codes)

# allocate patients to hospitals (uneven volumes)
hosp_weights <- rgamma(N_HOSPITALS, shape = 2, rate = 1); hosp_weights <- hosp_weights/sum(hosp_weights)
diag_hosp <- sample(hosp_codes, n, replace = TRUE, prob = hosp_weights)
diag_trust <- substr(diag_hosp, 1, 3)

# ----------------------------- 2. PATIENTS -----------------------------------
ydiag  <- as.integer(psamp(setNames(P$year_weights, P$year_levels), n))
# random diagnosis date within the calendar year
doy    <- as.integer(runif(n, 0, 364))
diagmdy<- as.Date(paste0(ydiag, "-01-01")) + doy

agediag <- round(pmin(99, pmax(20, rnorm(n, P$mean_age, P$sd_age))))
birthmdy<- diagmdy - round(agediag * 365.25)
sex     <- psamp(P$p_sex, n)
ethnicity_group_broad <- psamp(P$p_ethnicity, n)
imd     <- psamp(P$p_imd, n)
stage   <- psamp(P$p_stage, n)
cci_grp <- psamp(P$p_cci, n)
grade   <- psamp(P$p_grade, n)
route_combined <- psamp(P$p_route, n)
modality<- psamp(P$p_modality, n)
proc    <- psamp(P$p_proc, n)

# --------------------- 3. WAITING TIMES (multilevel) -------------------------
# Each component drawn from a Gamma with a mean shifted by hospital RE +
# covariate effects, so the downstream Stata models recover sensible signs and
# a non-trivial hospital-level ICC.
age_term <- P$age_slope * pmax(0, agediag - 70)
lp_extra <- P$eff_route[route_combined] + P$eff_cci[cci_grp] +
            P$eff_stage[stage] + P$year_effect[as.character(ydiag)] + age_term
lp_extra <- lp_extra - mean(lp_extra)   # centre: keep overall means on target, retain signal

draw_components <- function(idx) {
  md <- pmax(1, P$mean_dx_to_dtt + hosp_re_dtt[diag_hosp[idx]] + 0.45*lp_extra[idx])
  mt <- pmax(1, P$mean_dtt_to_tx + hosp_re_tx[diag_hosp[idx]]  + 0.55*lp_extra[idx])
  a  <- pmax(0L, as.integer(round(rgamma(length(idx), P$gamma_shape_dtt, scale = md/P$gamma_shape_dtt))))
  b  <- pmax(0L, as.integer(round(rgamma(length(idx), P$gamma_shape_tx,  scale = mt/P$gamma_shape_tx))))
  list(a = a, b = b)
}
wt_dx_to_dtt <- integer(n); wt_dtt_to_tx <- integer(n)
todo <- seq_len(n)
repeat {                                   # enforce total wait in [0,180] (cohort bound from 1b)
  cm <- draw_components(todo)
  wt_dx_to_dtt[todo] <- cm$a; wt_dtt_to_tx[todo] <- cm$b
  todo <- which((wt_dx_to_dtt + wt_dtt_to_tx) > 180)
  if (!length(todo)) break
}
wt_dx_to_tx  <- wt_dx_to_dtt + wt_dtt_to_tx          # identity preserved; total in [0,180]
days_diag_to_surg <- as.integer(wt_dx_to_tx)

# --------------------------- 4. PATHWAY DATES --------------------------------
dtt_date <- diagmdy + wt_dx_to_dtt
tx_date  <- diagmdy + wt_dx_to_tx
EPISTART <- tx_date
ADMIDATE <- tx_date
# MDT usually before DTT; ~15% missing
mdt_present <- runif(n) > 0.15
mdt_offset  <- round(runif(n) * wt_dx_to_dtt)
mdt_date    <- as.Date(ifelse(mdt_present, diagmdy + mdt_offset, NA), origin = "1970-01-01")
# CWT treatment date agrees with HES within +/-5 days
tx_jit      <- sample(-3:3, n, replace = TRUE)
cwt_tx_date <- tx_date + tx_jit
tx_date_diff<- as.numeric(abs(tx_date - cwt_tx_date))
# keep within the linkage window the real cohort enforces
fix <- tx_date_diff > 5
cwt_tx_date[fix] <- tx_date[fix]; tx_date_diff[fix] <- 0
first_hosp_date <- diagmdy - sample(0:21, n, replace = TRUE)

# ~10% of records have a missing decision-to-treat date in CWT. The real cohort
# only requires wt_dx_to_tx >= 0, so these keep the total wait but have NA for
# the DTT-based components (this exercises the Stata !missing() filters).
dtt_missing <- runif(n) < 0.10
dtt_date[dtt_missing]     <- NA
wt_dx_to_dtt[dtt_missing] <- NA
wt_dtt_to_tx[dtt_missing] <- NA
mdt_date[dtt_missing & !mdt_present] <- NA

dead <- as.integer(runif(n) < P$p_dead)
finmdy <- as.Date(ifelse(dead == 1, tx_date + sample(30:1600, n, replace = TRUE), NA),
                  origin = "1970-01-01")

# month/year helpers for CWT admin date fields
mon <- function(d) as.integer(format(d, "%m"))
yr  <- function(d) as.integer(format(d, "%Y"))
ds  <- function(d) format(d, "%d/%m/%Y")   # CWT stores some dates as dd/mm/YYYY strings

# --------------------- 5. PATHWAY ORDERING CHECK FLAGS -----------------------
dx_le_mdt_ok  <- is.na(mdt_date) | (diagmdy <= mdt_date)
dx_le_dtt_ok  <- is.na(dtt_date) | (diagmdy <= dtt_date)
dx_le_tx_ok   <-                    diagmdy <= tx_date
mdt_le_dtt_ok <- is.na(mdt_date) | is.na(dtt_date) | (mdt_date <= dtt_date)
dtt_le_tx_ok  <- is.na(dtt_date) |                    dtt_date <= tx_date
mdt_le_tx_ok  <- is.na(mdt_date) |                    mdt_date <= tx_date
seq_ok <- dx_le_dtt_ok & dx_le_mdt_ok & dx_le_tx_ok & mdt_le_dtt_ok & dtt_le_tx_ok & mdt_le_tx_ok

# ----------------------- 6. TUMOUR / CLINICAL DETAIL -------------------------
sub <- c(A="A", B="B", C="C")
stage_best <- paste0(stage, sample(c("","A","B","C"), n, TRUE, c(.5,.25,.18,.07)))
stage_best[stage_best %in% c("1B","1C")] <- "1A"   # tidy implausible substages
sitestr <- psamp(setNames(c(.05,.18,.17,.10,.20,.08,.07,.10,.05),
                          c("C18","C180","C181","C182","C183","C184","C185","C186","C189")), n)
site_icd10 <- sitestr                       # CWT site matches registry site -> site_match
site_match <- as.integer(substr(site_icd10,1,3) == substr(sitestr,1,3))
typestr <- psamp(setNames(c(.82,.08,.05,.05), c("8140","8480","8210","8000")), n)  # adenoca dominant
t_best <- psamp(setNames(c(.12,.30,.40,.18), paste0("T", 1:4)), n)
n_best <- psamp(setNames(c(.60,.28,.12),     paste0("N", 0:2)), n)
m_best <- rep("M0", n)                       # cohort excludes M1
dukes  <- ifelse(stage=="1","A", ifelse(stage=="2","B","C"))
nodesexcised <- as.integer(pmax(0, round(rnorm(n, 16, 7))))
nodesinvolved<- as.integer(pmin(nodesexcised, rpois(n, ifelse(stage=="3", 3, 0.4))))

# routes (raw fields consistent with route_combined)
final_route <- route_combined
route_bjc   <- route_combined
route_code  <- match(route_combined, names(P$p_route))
screendetected <- ifelse(route_combined == "Screening",
                         runif(n) < P$p_screen_given_route_screening,
                         runif(n) < 0.01)
tww_to_treat <- ifelse(route_combined == "Two Week Wait", wt_dx_to_tx, NA_integer_)

# OPCS for the primary resection
opcs_by_type <- list(
  right_hemi=c("H071","H072","H073"), left_hemi=c("H091","H092"),
  sigmoid=c("H101","H102","H103"), transverse=c("H081","H082"),
  total_subtotal=c("H041","H042","H051","H052"))
colon_opcs_primary <- vapply(proc, function(t) sample(opcs_by_type[[t]], 1), character(1))
opcs3 <- substr(colon_opcs_primary, 1, 3)

# treatment provider/site (mostly same trust as diagnosis; minority change trust)
change_trust <- runif(n) < 0.18
treat_trust  <- diag_trust
treat_trust[change_trust] <- sample(trust_codes, sum(change_trust), replace = TRUE)
SITETRET <- ifelse(change_trust,
                   paste0(treat_trust, sprintf("%02d", sample(1:9, n, replace=TRUE))),
                   diag_hosp)
PROCODE3 <- substr(SITETRET, 1, 3)
emergency <- runif(n) < P$p_emergency
ADMIMETH  <- ifelse(emergency, sample(c("21","22","24","2A"), n, replace=TRUE),
                    sample(c("11","12","13"), n, replace=TRUE))

# IDs and names
pid <- sprintf("P%08d", seq_len(n))

# ------------------- 7. ASSEMBLE ALL 172 COLUMNS (in order) ------------------
NA_chr <- rep(NA_character_, n); NA_num <- rep(NA_real_, n)
NA_date<- as.Date(rep(NA, n), origin = "1970-01-01")

cohort <- data.frame(
  pseudo_patientid = pid,
  pseudo_tumourid  = paste0(pid, "T1"),
  diagmdy          = diagmdy,
  ydiag            = ydiag,
  cancer           = factor("Colon"),
  sitestr          = sitestr,
  typestr          = typestr,
  basisofdiagnosis = psamp(setNames(c(.86,.07,.04,.03),
                            c("Histology","Cytology","Clinical","Investigation")), n),
  grade            = grade,
  stage_best       = stage_best,
  stage_best_system= "TNM8",
  t_best = t_best, n_best = n_best, m_best = m_best,
  t_path = t_best, n_path = n_best, m_path = m_best,
  sex              = factor(sex, levels = c("Male","Female")),
  agediag          = agediag,
  birthmdy         = birthmdy,
  ethnicity_group_broad = ethnicity_group_broad,
  lsoa11_code      = sprintf("E01%06d", sample(1:330000, n, replace = TRUE)),
  NHSE_reversed_imd_quintile_lsoas = as.integer(imd),
  canalliance_2024_code = sprintf("CA%02d", sample(1:21, n, replace = TRUE)),
  canalliance_2024_name = paste0("Cancer Alliance ", sample(1:21, n, replace = TRUE)),
  diag_trust       = diag_trust,
  diag_trust_name  = unname(trust_name[diag_trust]),
  first_trust      = diag_trust,
  first_trust_name = unname(trust_name[diag_trust]),
  first_hosp_date  = first_hosp_date,
  diag_hosp        = diag_hosp,
  diag_hosp_name   = unname(hosp_name[diag_hosp]),
  route_bjc        = route_bjc,
  final_route      = final_route,
  route_code       = route_code,
  tww_to_treat     = tww_to_treat,
  sg_flag          = 1L,
  rt_flag          = as.integer(runif(n) < 0.10),
  ct_flag          = as.integer(runif(n) < (ifelse(stage=="3",0.55,0.15))),
  screendetected   = as.integer(screendetected),
  dead             = dead,
  finmdy           = finmdy,
  dco              = 0L,
  er_status        = NA_chr, pr_status = NA_chr, her2_status = NA_chr,
  laterality.x     = NA_chr,
  dukes            = dukes,
  nodesexcised     = nodesexcised,
  nodesinvolved    = nodesinvolved,
  final_route_chr  = final_route,
  route_bjc_chr    = route_bjc,
  route_combined   = factor(route_combined, levels = names(P$p_route)),
  stage            = stage,
  STUDY_ID         = pid,
  ADMIDATE         = ADMIDATE,
  ADMIMETH         = ADMIMETH,
  PROCODE3         = PROCODE3,
  SITETRET         = SITETRET,
  EPISTART         = EPISTART,
  EPIORDER         = 1L,
  EPITYPE          = "1",
  stringsAsFactors = FALSE
)

# OPDATE_01..24 (Date): primary on OPDATE_01, rest NA
cohort$OPDATE_01 <- tx_date
for (i in 2:24) cohort[[sprintf("OPDATE_%02d", i)]] <- NA_date

cohort$op_position                 <- "OPERTN_01"
cohort$colon_opcs_primary          <- colon_opcs_primary
cohort$opcs3                       <- opcs3
cohort$colon_opcs_primary_position <- 1L
cohort$colon_proc_type_primary     <- proc
cohort$primary_flag                <- TRUE
cohort$all_colon_opcs              <- colon_opcs_primary
cohort$all_colon_proc_types        <- proc
cohort$all_colon_opcs_fields       <- "OPERTN_01"
cohort$n_colon_codes_in_episode    <- as.integer(1L + rbinom(n, 2, 0.08))
cohort$emergency                   <- emergency
cohort$days_diag_to_surg           <- days_diag_to_surg
cohort$tx_date                     <- tx_date

# CWT administrative block ----------------------------------------------------
cohort$org_ppi      <- diag_trust
cohort$ref_source   <- psamp(setNames(c(.5,.2,.15,.15),c("GP","Screening","Consultant","Other")), n)
cohort$priority_type<- psamp(setNames(c(.55,.3,.15),c("Urgent","Two Week Wait","Routine")), n)
cohort$dec_to_ref_date        <- diagmdy + sample(-3:3, n, replace = TRUE)
cohort$month_dec_to_ref_date  <- mon(cohort$dec_to_ref_date)
cohort$year_dec_to_ref_date   <- yr(cohort$dec_to_ref_date)
cohort$crtp_date              <- dtt_date
cohort$month_crtp_date        <- mon(dtt_date)
cohort$year_crtp_date         <- yr(dtt_date)
cohort$ref_type               <- psamp(setNames(c(.7,.3),c("Cancer","Non-cancer")), n)
cohort$cons_upgrade_date      <- NA_date
cohort$month_cons_upgrade_date<- NA_num
cohort$year_cons_upgrade_date <- NA_num
cohort$org_cons_upgrade       <- NA_chr
cohort$delay_cons_treat_reason<- NA_chr
cohort$date_first_seen        <- diagmdy + sample(0:14, n, replace = TRUE)
cohort$month_date_first_seen  <- mon(cohort$date_first_seen)
cohort$year_date_first_seen   <- yr(cohort$date_first_seen)
cohort$org_first_seen         <- diag_trust
cohort$wta_first_seen         <- as.integer(pmax(0, round(rnorm(n, 7, 5))))
cohort$wta_first_seen_reason  <- NA_chr
cohort$delay_ref_fs_reason    <- NA_chr
cohort$patient_status         <- psamp(setNames(c(.95,.05),c("NHS","Private")), n)
cohort$laterality.y           <- NA_chr
cohort$mets_site              <- NA_chr     # cohort is M0
cohort$treat_period_start     <- ds(dtt_date)     # dd/mm/YYYY string (as in CWT)
cohort$month_treat_period_start <- mon(dtt_date)
cohort$year_treat_period_start  <- yr(dtt_date)
cohort$org_dec_to_treat       <- treat_trust
cohort$fdp_end_reason         <- psamp(setNames(c(.9,.1),c("First treatment","Other")), n)
cohort$fdp_diag_site          <- diag_hosp
cohort$fdp_end_date           <- tx_date
cohort$month_fdp_end_date     <- mon(tx_date)
cohort$year_fdp_end_date      <- yr(tx_date)
cohort$delay_fdp_reason       <- NA_chr
cohort$fdp_exclusion_reason   <- NA_chr
cohort$fdp_outcome_prof_type  <- psamp(setNames(c(.8,.2),c("Surgeon","Other")), n)
cohort$fdp_outcome_method     <- "Surgery"
cohort$org_fdp_end            <- treat_trust
cohort$treat_start            <- ds(cwt_tx_date)  # dd/mm/YYYY string (as in CWT)
cohort$month_treat_start      <- mon(cwt_tx_date)
cohort$year_treat_start       <- yr(cwt_tx_date)
cohort$org_treat_start        <- SITETRET         # used as trust_treat fallback in Stata
cohort$cte_type               <- "Surgery"
cohort$modality               <- modality
cohort$clin_trial             <- as.integer(runif(n) < 0.03)
cohort$care_setting           <- psamp(setNames(c(.97,.03),c("NHS","Private")), n)
cohort$delay_dtt_treat_reason <- NA_chr
cohort$wta_treat              <- as.integer(pmax(0, round(rnorm(n, 14, 8))))
cohort$wta_treat_reason       <- NA_chr
cohort$delay_ref_treat_reason <- NA_chr
cohort$radio_intent           <- NA_chr
cohort$radio_priority         <- NA_chr
cohort$mdt_ind                <- ifelse(mdt_present, "Y", "N")
cohort$mdt_date               <- mdt_date
cohort$month_mdt_date         <- mon(mdt_date)
cohort$year_mdt_date          <- yr(mdt_date)
cohort$practice_code          <- sprintf("%s%05d", sample(LETTERS, n, replace=TRUE),
                                         sample(10000:99999, n, replace=TRUE))
cohort$site_icd10             <- site_icd10
cohort$dtt_date               <- dtt_date
cohort$cwt_tx_date            <- cwt_tx_date
cohort$tx_date_diff           <- tx_date_diff
cohort$dx_le_mdt_ok  <- dx_le_mdt_ok
cohort$dx_le_dtt_ok  <- dx_le_dtt_ok
cohort$dx_le_tx_ok   <- dx_le_tx_ok
cohort$mdt_le_dtt_ok <- mdt_le_dtt_ok
cohort$dtt_le_tx_ok  <- dtt_le_tx_ok
cohort$mdt_le_tx_ok  <- mdt_le_tx_ok
cohort$seq_ok        <- seq_ok
cohort$wt_dx_to_dtt  <- as.numeric(wt_dx_to_dtt)
cohort$wt_dx_to_tx   <- as.numeric(wt_dx_to_tx)
cohort$wt_dtt_to_tx  <- as.numeric(wt_dtt_to_tx)
cohort$site_match    <- site_match

# --------------------------- 8. STRUCTURE CHECK ------------------------------
target_names <- c(
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

cohort <- cohort[, target_names]   # enforce exact order
stopifnot(ncol(cohort) == 172L)
stopifnot(identical(names(cohort), target_names))
ok_id <- with(cohort, wt_dx_to_tx == wt_dx_to_dtt + wt_dtt_to_tx)
stopifnot(all(ok_id[!is.na(ok_id)]))   # identity holds wherever DTT is observed

cat(sprintf("Synthetic cohort built: %d rows x %d cols across %d hospitals / %d trusts\n",
            nrow(cohort), ncol(cohort), N_HOSPITALS, n_trust))

# ------------------------------ 9. WRITE OUT ---------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
rds_path <- file.path(OUT_DIR, "synthetic_colon_cohort.rds")
dta_path <- file.path(OUT_DIR, "synthetic_colon_cohort.dta")
saveRDS(cohort, rds_path)

write_dta_safe <- function(df, path) {
  # Stata variable names cannot contain '.', so laterality.x/.y -> laterality_x/_y
  # (this is the same conversion Stata applies on import). The .rds keeps the
  # exact original R names; only the .dta is sanitised.
  d <- df
  names(d) <- gsub("[^A-Za-z0-9_]", "_", names(d))
  if (requireNamespace("haven", quietly = TRUE)) {
    haven::write_dta(d, path); return("haven")
  }
  # foreign fallback: coerce types Stata-friendly
  for (j in seq_along(d)) {
    x <- d[[j]]
    if (inherits(x, "Date")) d[[j]] <- as.numeric(x) - as.numeric(as.Date("1960-01-01")) # Stata epoch
    else if (is.logical(x))  d[[j]] <- as.integer(x)
    else if (is.factor(x))   d[[j]] <- as.character(x)
  }
  foreign::write.dta(d, path, version = 10L); return("foreign")
}
eng <- write_dta_safe(cohort, dta_path)
cat("Wrote:", rds_path, "\n      ", dta_path, sprintf("(via %s)\n", eng))
