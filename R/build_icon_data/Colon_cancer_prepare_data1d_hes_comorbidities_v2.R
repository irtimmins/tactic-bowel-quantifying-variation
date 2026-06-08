


library(tidyverse)

colon_cohort <-   readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_2015_2022.rds")


hes_apc_raw <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/hes_apc_raw_colon_2014_2022.rds")

# A. RCS CHARLSON ICD-10 CODES (Armitage & van der Meulen 2010, Table 1) 
# Using 3-digit prefixes as per RCS coding philosophy
# Dots removed to match DIAG_4_ format in HES

# Conditions to look for in PREVIOUS admissions only (and index admission
# for non-acute codes). Acute conditions (*) from previous admissions only.

charlson_icd10 <- list(
  myocardial_infarction    = list(
    codes = c("I21","I22","I23","I252"),
    acute = TRUE    # I21, I22, I23 are acute - previous admissions only
  ),
  congestive_heart_failure = list(
    codes = c("I11","I13","I255","I42","I43","I50","I517"),
    acute = FALSE
  ),
  peripheral_vascular      = list(
    codes = c("I70","I71","I72","I73","I770","I771","K551","K558","K559","R02","Z958","Z959"),
    acute = FALSE
  ),
  cerebrovascular          = list(
    codes = c("G45","G46","I60","I61","I62","I63","I64","I65","I66","I67","I68","I69"),
    acute = FALSE
  ),
  dementia                 = list(
    codes = c("A810","F00","F01","F02","F03","F051","G30","G31"),
    acute = FALSE
  ),
  chronic_pulmonary        = list(
    codes = c("I26","I27","J40","J41","J42","J43","J44","J45","J47",
              "J60","J61","J62","J63","J64","J65","J66","J67","J684","J701","J703"),
    acute = FALSE
    # J46 (status asthmaticus) is acute - omitted as previous-only and rare
  ),
  rheumatological          = list(
    codes = c("M05","M06","M09","M120","M315","M32","M33","M34","M35","M36"),
    acute = FALSE
  ),
  liver                    = list(
    codes = c("B18","I85","I864","I982","K70","K71","K721","K729","K76","R162","Z944"),
    acute = FALSE
  ),
  diabetes                 = list(
    codes = c("E10","E11","E12","E13","E14"),
    acute = FALSE
  ),
  hemiplegia_paraplegia    = list(
    codes = c("G114","G81","G82","G83"),
    acute = FALSE
  ),
  renal                    = list(
    codes = c("I12","I13","N01","N03","N05","N07","N08","N18","N19","N25","Z49","Z940","Z992"),
    acute = FALSE
    # N171, N172 acute - previous admissions only, already handled
  ),
  malignancy               = list(
    codes = c("C00","C01","C02","C03","C04","C05","C06","C07","C08","C09",
              "C10","C11","C12","C13","C14","C15","C16","C17","C18","C19",
              "C20","C21","C22","C23","C24","C25","C26",
              "C30","C31","C32","C33","C34","C37","C38","C39","C40","C41",
              "C43","C45","C46","C47","C48","C49","C50","C51","C52","C53",
              "C54","C55","C56","C57","C58",
              "C60","C61","C62","C63","C64","C65","C66","C67","C68","C69",
              "C70","C71","C72","C73","C74","C75","C76",
              "C80","C81","C82","C83","C84","C85","C88","C90","C91","C92",
              "C93","C94","C95","C96","C97"),
    acute = FALSE
  ),
  metastatic_solid_tumour  = list(
    codes = c("C77","C78","C79"),
    acute = FALSE
  ),
  aids                     = list(
    codes = c("B20","B21","B22","B24"),
    acute = FALSE
  )
)

# Flat lookup for joining
charlson_lookup <- imap_dfr(charlson_icd10, function(x, condition) {
  tibble(prefix = x$codes, condition = condition, acute = x$acute)
})



# Diagnosis fields to use (positions 1-7 only, per RCS paper)
# DIAG_4_01 = primary diagnosis included throughout
sec_diag_cols <- paste0("DIAG_4_", str_pad(1:20, 2, pad = "0"))
sec_diag_cols 

# B. INDEX ADMISSION EPISODES 
# Non-acute conditions from the index admission episode
# Index admission = the surgical episode in colon_cohort (tx_date or ADMIDATE)

index_episodes <- colon_cohort %>%
  select(STUDY_ID = pseudo_patientid, tx_date) %>%
  mutate(STUDY_ID = as.character(STUDY_ID),
         tx_date  = as.Date(tx_date)) %>%
  inner_join(
    hes_apc_raw %>% select(STUDY_ID, ADMIDATE, any_of(sec_diag_cols)),
    by = "STUDY_ID"
  ) %>%
  filter(ADMIDATE == tx_date)   # match on admission date of surgical episode

# C. PREVIOUS ADMISSIONS (12-month lookback) 

diag_dates <- colon_cohort %>%
  select(STUDY_ID = pseudo_patientid, diagmdy) %>%
  mutate(
    STUDY_ID       = as.character(STUDY_ID),
    diagmdy        = as.Date(diagmdy),
    lookback_start = diagmdy - 365,
    lookback_end   = diagmdy - 1
  )

prev_episodes <- hes_apc_raw %>%
  select(STUDY_ID, EPISTART, any_of(sec_diag_cols)) %>%
  inner_join(diag_dates, by = "STUDY_ID") %>%
  filter(EPISTART >= lookback_start, EPISTART <= lookback_end) %>%
  select(STUDY_ID, any_of(sec_diag_cols))





# How many patients have at least one HES admission in the 12 months before diagnosis?
hes_apc_raw %>%
  inner_join(diag_dates, by = "STUDY_ID") %>%
  filter(EPISTART >= lookback_start, EPISTART <= lookback_end) %>%
  distinct(STUDY_ID) %>%
  nrow()

# D. PIVOT AND MATCH 

pivot_and_match <- function(episodes, source_label) {
  episodes %>%
    pivot_longer(
      cols      = any_of(sec_diag_cols),
      names_to  = "diag_position",
      values_to = "icd_code"
    ) %>%
    filter(!is.na(icd_code), icd_code != "-", icd_code != "") %>%
    mutate(
      icd_code = str_remove_all(str_trim(icd_code), "\\."),
      source   = source_label,
      p3 = str_sub(icd_code, 1, 3),
      p4 = str_sub(icd_code, 1, 4)
    ) %>%
    left_join(charlson_lookup %>% filter(nchar(prefix) == 3) %>%
                rename(cond3 = condition, acute3 = acute),
              by = c("p3" = "prefix"),
              relationship = "many-to-many") %>%
    left_join(charlson_lookup %>% filter(nchar(prefix) == 4) %>%
                rename(cond4 = condition, acute4 = acute),
              by = c("p4" = "prefix"),
              relationship = "many-to-many") %>%
    mutate(
      charlson_condition = coalesce(cond4, cond3),
      is_acute           = coalesce(acute4, acute3)
    ) %>%
    filter(!is.na(charlson_condition)) %>%
    select(STUDY_ID, charlson_condition, is_acute, source)
}

prev_matched  <- pivot_and_match(prev_episodes,  "previous")
index_matched <- pivot_and_match(index_episodes, "index") %>%
  filter(!is_acute)   # acute conditions excluded from index admission

#all_matched <- bind_rows(prev_matched, index_matched)
all_matched <- prev_matched # Depends on methods whether to include index.

cat("Charlson-relevant rows (previous):", nrow(prev_matched), "\n")
cat("Charlson-relevant rows (index):   ", nrow(index_matched), "\n")

# E. COLLAPSE TO PATIENT LEVEL 

cci_patient <- all_matched %>%
  distinct(STUDY_ID, charlson_condition) %>%
  group_by(STUDY_ID) %>%
  summarise(
    cci_n_conditions = n_distinct(charlson_condition),
    cci_conditions   = paste(sort(unique(charlson_condition)), collapse = "; "),
    .groups = "drop"
  )

all_patients <- tibble(STUDY_ID = as.character(unique(colon_cohort$pseudo_patientid)))

cci_patient <- all_patients %>%
  left_join(cci_patient, by = "STUDY_ID") %>%
  mutate(
    cci_n_conditions = replace_na(cci_n_conditions, 0),
    cci_conditions   = replace_na(cci_conditions, "none"),
    cci_group = factor(
      case_when(
        cci_n_conditions == 0 ~ "0",
        cci_n_conditions == 1 ~ "1",
        cci_n_conditions == 2 ~ "2",
        TRUE                  ~ "3+"
      ),
      levels = c("0", "1", "2", "3+")
    )
  )

cat("\nRCS Charlson CCI distribution:\n")
print(table(cci_patient$cci_group))

colon_cohort <- colon_cohort %>%
  left_join(cci_patient, by = c("pseudo_patientid" = "STUDY_ID"))


colon_cohort %>%
  mutate(age_group = cut(agediag, breaks = c(0,49,59,69,79,Inf),
                         labels = c("<50","50-59","60-69","70-79","80+"))) %>%
  count(age_group, cci_group) %>%
  group_by(age_group) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  select(-n) %>%
  pivot_wider(names_from = cci_group, values_from = pct)

all_matched %>%
  filter(charlson_condition == "malignancy") %>%
  distinct(STUDY_ID) %>%
  nrow()


names(colon_cohort)


saveRDS(
  colon_cohort,
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_cci_2015_2022.rds"
)



