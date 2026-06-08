
# ---------------------------------------------------------------------
# 6. Build NCRAS anchor: one row per patient, diagnosis date and covariates
# ---------------------------------------------------------------------

ncras_index <- ncras_colon %>%
  mutate(
    pseudo_patientid = as.character(pseudo_patientid),
    diagmdy          = as.Date(diagmdy)
  ) %>%
  distinct(pseudo_patientid, .keep_all = TRUE)

# ---------------------------------------------------------------------
# 7. Join HES index surgery to NCRAS anchor
#    - elective episodes only
#    - plausible timing window: -90 to +180 days from diagnosis
#    - one row per patient: earliest post-dx elective resection
# ---------------------------------------------------------------------

ncras_hes <- ncras_index %>%
  left_join(
    hes_colon_episodes %>%
#      filter(!emergency) %>%
      mutate(pseudo_patientid = STUDY_ID),
    by = "pseudo_patientid"
  ) %>%
  mutate(days_diag_to_surg = as.integer(EPISTART - diagmdy)) %>%
  filter(!is.na(days_diag_to_surg),
         days_diag_to_surg >= -90,
         days_diag_to_surg <= 180) %>%
  # One row per patient: earliest plausible post-dx elective resection
  arrange(pseudo_patientid, EPISTART, EPIORDER) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  mutate(tx_date = EPISTART)

# QA
cat("NCRAS patients:                  ", n_distinct(ncras_index$pseudo_patientid), "\n")
cat("NCRAS patients with HES surgery: ", n_distinct(ncras_hes$pseudo_patientid), "\n")
cat("Missing HES match:               ",
    sum(is.na(ncras_hes$EPISTART)), "\n")

# Waiting time summary after HES join
ncras_hes %>%
  filter(!is.na(days_diag_to_surg)) %>%
  summarise(
    n           = n(),
    p25         = quantile(days_diag_to_surg, 0.25),
    mean        = mean(days_diag_to_surg),
    median      = median(days_diag_to_surg),
    p75         = quantile(days_diag_to_surg, 0.75),
    p90         = quantile(days_diag_to_surg, 0.90),
    max         = max(days_diag_to_surg),
    pct_over_62 = round(100 * mean(days_diag_to_surg > 62), 1),
    pct_neg     = round(100 * mean(days_diag_to_surg < 0), 1),
    pct_zero    = round(100 * mean(days_diag_to_surg == 0), 1)
  ) %>% print()

sum(ncras_hes$days_diag_to_surg == 0)

# ---------------------------------------------------------------------
# 8. Join CWT to NCRAS-HES and build final analysis cohort
#    - match on patient ID and treatment date (within 5 days)
#    - validate pathway ordering: Dx <= MDT <= DTT <= Tx
#    - deduplicate: prefer closest HES-CWT date match, then earliest DTT
# ---------------------------------------------------------------------

cwt_colon <- open_dataset(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/11_CWT_data_partitioned"
) %>%
  filter(site_icd10 %in% colon_icd) %>%
  collect()

#summary(as.factor(cwt_colon$site_icd10))

names(cwt_colon)

library(lubridate)

#View(cwt_colon)

cwt_colon %>%
  # filter(modality %in% c(1, 23, 24)) %>%
  mutate(
    dtt_date    = as.Date(treat_period_start, format = "%d/%m/%Y"),  # <-- DTT confirmed
    cwt_tx_date = as.Date(treat_start,        format = "%d/%m/%Y"),
    mdt_date    = as.Date(mdt_date,           format = "%d/%m/%Y")
  ) %>%
  pull(  dtt_date) %>% summary()

summary(cwt_colon$treat_period_start)

cwt_colon %>%
  mutate(year = year(as.Date(treat_start, format = "%d/%m/%Y"))) %>%
  count(year, modality) %>%
  arrange(year, modality) %>%
  print(n = Inf)

cwt_colon %>%
  mutate(year = year(as.Date(treat_start, format = "%d/%m/%Y"))) %>%
  count(year, modality) %>%
  group_by(year) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(year, modality) %>%
  print(n = Inf)


cwt_colon_surgery <-  cwt_colon %>%
  filter(
    modality %in% c("01", "23", "24"),
    # Exclude early 2020 modality 23/24 records prior to the coding transition
    !(modality %in% c("23", "24") & 
        as.Date(treat_start, format = "%d/%m/%Y") < as.Date("2020-06-01"))
  ) %>% 
  mutate(
    dtt_date    = as.Date(treat_period_start, format = "%d/%m/%Y"),  # <-- DTT confirmed
    cwt_tx_date = as.Date(treat_start,        format = "%d/%m/%Y"),
    mdt_date    = as.Date(mdt_date,           format = "%d/%m/%Y")
  )




colon_cohort <- ncras_hes %>%
  filter(!is.na(tx_date)) %>%
  left_join(cwt_colon_surgery,
    by = "pseudo_patientid"
  ) %>%
  # Compare HES-APC treatment date with CWT to confirm linkage
  mutate(tx_date_diff = as.numeric(abs(tx_date - cwt_tx_date))) %>%
  filter(!is.na(tx_date_diff), tx_date_diff <= 5) %>%
#  filter(!is.na(tx_date_diff), tx_date_diff == 0) %>%
  mutate(
    # Pairwise ordering checks; NA = can't check, not a failure
    dx_le_mdt_ok  = is.na(mdt_date) | (diagmdy <= mdt_date),
    dx_le_dtt_ok  = is.na(dtt_date) | (diagmdy <= dtt_date),
    dx_le_tx_ok   =                    diagmdy <= tx_date,
    mdt_le_dtt_ok = is.na(mdt_date) | is.na(dtt_date) | (mdt_date <= dtt_date),
    dtt_le_tx_ok  = is.na(dtt_date) |                    dtt_date <= tx_date,
    mdt_le_tx_ok  = is.na(mdt_date) |                    mdt_date <= tx_date,
    # Overall plausibility: Dx <= MDT <= DTT <= Tx (where observed)
    seq_ok = dx_le_dtt_ok & dx_le_mdt_ok & dx_le_tx_ok
    & mdt_le_dtt_ok & dtt_le_tx_ok & mdt_le_tx_ok
  ) %>%
  mutate(
    # Waiting time intervals (days)
    wt_dx_to_dtt  = as.numeric(dtt_date - diagmdy),
    # wt_dx_to_mdt  = as.numeric(mdt_date - diagmdy),
    wt_dx_to_tx   = as.numeric(tx_date  - diagmdy),
   # wt_mdt_to_dtt = as.numeric(dtt_date - mdt_date),
    wt_dtt_to_tx  = as.numeric(tx_date  - dtt_date)
  ) %>%
     # filter(
     #   if_all(starts_with("wt_"),
     #         ~ is.na(.x) | (.x >= 0 & .x <= 180))
     # ) %>%
  filter(wt_dx_to_tx >= 0 ) %>%
#  filter(wt_dx_to_dtt >= 0, wt_dtt_to_tx >= 0) %>%
 
   # Deduplicate: prefer closest HES-CWT date match, then earliest DTT
  mutate(site_match = as.integer(
    str_sub(site_icd10, 1, 3) == str_sub(sitestr, 1, 3)
  )) %>% 
  arrange(pseudo_patientid, desc(site_match), tx_date_diff, dtt_date, cwt_tx_date) %>%
 # group_by(pseudo_patientid) %>%
  #filter(n() > 1) %>%
  #View()
  distinct(pseudo_patientid, .keep_all = TRUE)

length(unique(colon_cohort$pseudo_patientid))  
  
summary(as.factor(colon_cohort$ydiag))

sum(colon_cohort$days_diag_to_surg == 0)
sum(colon_cohort$wt_dx_to_dtt == 0)
sum(colon_cohort$wt_dx_to_tx == 0)

summary(colon_cohort$wt_dx_to_tx)
hist(colon_cohort$wt_dx_to_tx, xlim = c(0, 180), breaks = 100)
sum(colon_cohort$wt_dx_to_tx == 0, na.rm = T)
summary(colon_cohort$wt_dx_to_tx)
summary(colon_cohort$wt_dx_to_dtt)

sum(colon_cohort$tx_date_diff == 0, na.rm = T)/length(colon_cohort$tx_date_diff)
sum(!is.na(colon_cohort$tx_date_diff))
sum(is.na(colon_cohort$tx_date_diff))
summary(colon_cohort$tx_date)
summary(colon_cohort$cwt_tx_date)
hist(colon_cohort$tx_date_diff, breaks = 10000, xlim = c(-5, 50))
hist(colon_cohort$wt_dx_to_dtt, breaks = 10000, xlim = c(0, 180))


summary(colon_cohort$wt_dtt_to_tx)
sd(colon_cohort$wt_dtt_to_tx[colon_cohort$wt_dtt_to_tx>=0])


summary(colon_cohort$wt_dx_to_dtt[colon_cohort$wt_dx_to_dtt >= 0])
sd(colon_cohort$wt_dx_to_dtt[colon_cohort$wt_dx_to_dtt >= 0])

summary(colon_cohort$wt_dx_to_dtt[colon_cohort$wt_dx_to_dtt >= 0])

summary(colon_cohort$wt_dx_to_tx)

summary(colon_cohort$wt_dtt_to_tx)

sd(colon_cohort$wt_dtt_to_tx)
#8288/54749

saveRDS(
  colon_cohort,
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_2015_2022.rds"
)



summary(colon_cohort$tx_date_diff)
hist(colon_cohort$tx_date_diff, breaks = 1000, xlim = c(0,20))
sum(colon_cohort$tx_date_diff == 0, na.rm = T)
sum(abs(colon_cohort$tx_date_diff) < 5, na.rm = T)
summary(colon_cohort$tx_date_diff)
summary(colon_cohort$cwt_tx_date)

test <- colon_cohort %>%
  filter(nchar(diag_hosp) == 5) %>%
  group_by(diag_hosp) %>%
  summarise(count = n(),
            est = mean(wt_dx_to_dtt)) %>%
  filter(count > 50)

summary(test$count)
#summary(as.factor(test$diag_hosp))
summary(test$est)
mean(test$est)
mean(colon_cohort$wt_dx_to_dtt)

length(unique(colon_cohort$diag_hosp))
colon_cohort$diag_hosp

colon_cohort_qc <- colon_cohort %>%
  mutate(diag_hosp = trimws(as.character(diag_hosp))) %>%
  filter(grepl("^R[A-Z0-9]{4}$", diag_hosp))

length(unique(colon_cohort_qc$diag_hosp))
sort(unique(colon_cohort$diag_hosp))


sum(colon_cohort$diag_hosp == "")

colon_cohort_qc
sum(grepl("^R[A-Z0-9]{4}$", colon_cohort$diag_hosp))
