# ---------------------------------------------------------------------
# CHARLSON COMORBIDITY INDEX (NBOCA method)
# Quan et al. (2005) ICD-10 coding algorithms
# Secondary diagnoses (DIAG_4_02 to DIAG_4_20) only
# 12 months prior to diagnosis date
# Grouped: 0 / 1 / 2+
# ---------------------------------------------------------------------

colon_cohort <-   readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_2015_2022.rds")

# A. ICD-10 CODE DEFINITIONS (Quan 2005) 

charlson_icd10 <- list(
  myocardial_infarction    = c("I21","I22","I252"),
  congestive_heart_failure = c("I099","I110","I130","I132","I255",
                               "I420","I425","I426","I427","I428",
                               "I429","I43","I50","P290"),
  peripheral_vascular      = c("I70","I71","I731","I738","I739",
                               "I771","I790","I792","K551","K558",
                               "K559","Z958","Z959"),
  cerebrovascular          = c("G45","G46","H340","I60","I61","I62",
                               "I63","I64","I65","I66","I67","I68","I69"),
  dementia                 = c("F00","F01","F02","F03","F09","G30",
                               "G311","G312","G315","G318","G319",
                               "G328","G914","G948","R413","R54"),
  chronic_pulmonary        = c("I278","I279","J40","J41","J42","J43",
                               "J44","J45","J46","J47","J60","J61","J62",
                               "J63","J64","J65","J66","J67","J684",
                               "J701","J703"),
  rheumatic                = c("M05","M06","M315","M32","M33","M34",
                               "M351","M353","M360"),
  peptic_ulcer             = c("K25","K26","K27","K28"),
  mild_liver               = c("B18","K700","K701","K702","K703",
                               "K709","K713","K714","K715","K717",
                               "K73","K74","K760","K762","K763",
                               "K764","K768","K769","Z944"),
  diabetes_uncomplicated   = c("E100","E101","E106","E108","E109",
                               "E110","E111","E116","E118","E119",
                               "E120","E121","E126","E128","E129",
                               "E130","E131","E136","E138","E139",
                               "E140","E141","E146","E148","E149"),
  diabetes_complicated     = c("E102","E103","E104","E105","E107",
                               "E112","E113","E114","E115","E117",
                               "E122","E123","E124","E125","E127",
                               "E132","E133","E134","E135","E137",
                               "E142","E143","E144","E145","E147"),
  hemiplegia               = c("G041","G114","G801","G802","G81",
                               "G82","G830","G831","G832","G833",
                               "G834","G839"),
  renal                    = c("I120","I131","N032","N033","N034",
                               "N035","N036","N037","N052","N053",
                               "N054","N055","N056","N057","N18",
                               "N19","N250","Z490","Z491","Z492",
                               "Z940","Z992"),
  severe_liver             = c("I850","I859","I864","I982","K704",
                               "K711","K721","K729","K765","K766",
                               "K767"),
  aids                     = c("B20","B21","B22","B24")
)

# Build flat lookup: ICD prefix -> condition name
charlson_lookup <- imap_dfr(charlson_icd10, function(codes, condition) {
  tibble(prefix = codes, condition = condition)
})


#  B. DIAGNOSIS DATE LOOKBACK WINDOW 

diag_dates <- colon_cohort %>%
  select(STUDY_ID = pseudo_patientid, diagmdy) %>%
  mutate(
    STUDY_ID       = as.character(STUDY_ID),
    diagmdy        = as.Date(diagmdy),
    lookback_start = diagmdy - 365,
    lookback_end   = diagmdy - 1
  )


# C. FILTER HES TO SECONDARY DIAGNOSES IN LOOKBACK WINDOW 
# DIAG_4_01 = primary diagnosis - exclude it, use _02 to _20 only

sec_diag_cols <- paste0("DIAG_4_", str_pad(2:20, 2, pad = "0"))

hes_apc_raw <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/hes_apc_raw_colon_2014_2022.rds")

hes_lookback <- hes_apc_raw %>%
  select(STUDY_ID, EPISTART, ADMIDATE, any_of(sec_diag_cols)) %>%
  inner_join(diag_dates, by = "STUDY_ID") %>%
  filter(EPISTART >= lookback_start, EPISTART <= lookback_end) %>%
  select(STUDY_ID, EPISTART, any_of(sec_diag_cols))

# Spot check - what do the raw secondary diag codes look like?
hes_lookback %>% 
  select(STUDY_ID, EPISTART, DIAG_4_02, DIAG_4_03) %>% 
  filter(!is.na(DIAG_4_02)) %>% 
  head(20)

# D. PIVOT LONG, CLEAN, MATCH TO CHARLSON CONDITIONS 

hes_diag_long <- hes_lookback %>%
  pivot_longer(
    cols      = any_of(sec_diag_cols),
    names_to  = "diag_position",
    values_to = "icd_code"
  ) %>%
  filter(!is.na(icd_code), icd_code != "-", icd_code != "") %>%
  mutate(
    # DIAG_4_ cols contain codes like "I500" - strip any dots just in case
    icd_code = str_remove_all(str_trim(icd_code), "\\.")
  )

# Vectorised prefix match via a join on first N characters
# Much faster than row-wise map_chr for large datasets
match_charlson_join <- function(diag_long, lookup) {
  # Try matching on 3-char prefix first, then 4-char
  # Build all possible prefixes from each code
  diag_long %>%
    mutate(
      p3 = str_sub(icd_code, 1, 3),
      p4 = str_sub(icd_code, 1, 4)
    ) %>%
    left_join(lookup %>% filter(nchar(prefix) == 3) %>% rename(cond3 = condition),
              by = c("p3" = "prefix")) %>%
    left_join(lookup %>% filter(nchar(prefix) == 4) %>% rename(cond4 = condition),
              by = c("p4" = "prefix")) %>%
    mutate(
      # 4-char match takes priority (more specific)
      charlson_condition = coalesce(cond4, cond3)
    ) %>%
    select(-p3, -p4, -cond3, -cond4)
}

hes_diag_long <- match_charlson_join(hes_diag_long, charlson_lookup) %>%
  filter(!is.na(charlson_condition))

cat("Charlson-relevant diagnosis rows:", nrow(hes_diag_long), "\n")

# What codes are actually being matched?
head(sort(table(hes_diag_long$charlson_condition), decreasing = TRUE), 20)


# E. COLLAPSE TO PATIENT LEVEL 

cci_patient <- hes_diag_long %>%
  distinct(STUDY_ID, charlson_condition) %>%
  group_by(STUDY_ID) %>%
  summarise(
    cci_n_conditions = n_distinct(charlson_condition),
    cci_conditions   = paste(sort(unique(charlson_condition)), collapse = "; "),
    .groups = "drop"
  )

# Left join to all patients so those with no comorbidities get score 0
cci_patient <- tibble(STUDY_ID = as.character(unique(colon_cohort$pseudo_patientid))) %>%
  left_join(cci_patient, by = "STUDY_ID") %>%
  mutate(
    cci_n_conditions = replace_na(cci_n_conditions, 0),
    cci_conditions   = replace_na(cci_conditions, "none"),
    cci_group        = factor(
      case_when(
        cci_n_conditions == 0 ~ "0",
        cci_n_conditions == 1 ~ "1",
        TRUE                  ~ "2+"
      ),
      levels = c("0", "1", "2+")
    )
  )

cat("\nCharlson CCI distribution:\n")
print(table(cci_patient$cci_group))


# F. MERGE BACK TO COHORT  

colon_cohort <- colon_cohort %>%
  left_join(cci_patient, by = c("pseudo_patientid" = "STUDY_ID"))

cat("\nCCI merged. Check:\n")
print(table(colon_cohort$cci_group, useNA = "ifany"))
