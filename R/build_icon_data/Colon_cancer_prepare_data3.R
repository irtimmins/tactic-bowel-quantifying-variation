

# =============================================================================
# EXPLORATORY ANALYSIS - Patient and tumour characteristics
# =============================================================================

colon_cohort <- readRDS(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_2015_2022.rds"
)


# -----------------------------------------------------------------------------
# 1. Age at diagnosis
# -----------------------------------------------------------------------------

cat("=== AGE AT DIAGNOSIS ===\n")
colon_cohort %>%
  summarise(
    n      = n(),
    mean   = round(mean(agediag,   na.rm = TRUE), 1),
    sd     = round(sd(agediag,     na.rm = TRUE), 1),
    median = median(agediag,       na.rm = TRUE),
    p25    = quantile(agediag, 0.25, na.rm = TRUE),
    p75    = quantile(agediag, 0.75, na.rm = TRUE),
    min    = min(agediag,          na.rm = TRUE),
    max    = max(agediag,          na.rm = TRUE)
  ) %>% print()

# Age groups
colon_cohort %>%
  mutate(age_group = cut(agediag,
                         breaks = c(0, 49, 59, 69, 79, Inf),
                         labels = c("<50", "50-59", "60-69", "70-79", "80+"),
                         right  = TRUE)) %>%
  count(age_group) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

# Age by year - detect any drift
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  group_by(year_diag) %>%
  summarise(
    n           = n(),
    median_age  = median(agediag, na.rm = TRUE),
    mean_age    = round(mean(agediag, na.rm = TRUE), 1),
    pct_over_80 = round(100 * mean(agediag > 80, na.rm = TRUE), 1)
  ) %>%
  print()





# -----------------------------------------------------------------------------
# 2. Sex
# -----------------------------------------------------------------------------

cat("\n=== SEX ===\n")
colon_cohort %>%
  count(sex) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

# Sex by year
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  count(year_diag, sex) %>%
  group_by(year_diag) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(year_diag, sex) %>%
  print()

# Waiting times by sex
colon_cohort %>%
  group_by(sex) %>%
  summarise(
    n             = n(),
    median_dx_dtt = median(wt_dx_to_dtt, na.rm = TRUE),
    median_dtt_tx = median(wt_dtt_to_tx, na.rm = TRUE),
    median_dx_tx  = median(wt_dx_to_tx,  na.rm = TRUE),
    pct_over_62   = round(100 * mean(wt_dx_to_tx > 62,  na.rm = TRUE), 1),
    pct_over_104  = round(100 * mean(wt_dx_to_tx > 104, na.rm = TRUE), 1)
  ) %>%
  print()

# -----------------------------------------------------------------------------
# 3. Ethnicity
# -----------------------------------------------------------------------------

cat("\n=== ETHNICITY ===\n")
colon_cohort %>%
  count(ethnicity_group_broad) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n)) %>%
  print()

# Missing ethnicity by year
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  group_by(year_diag) %>%
  summarise(
    n              = n(),
    pct_miss_ethn  = round(100 * mean(is.na(ethnicity_group_broad) |
                                        ethnicity_group_broad == "Unknown"), 1)
  ) %>%
  print()

# Waiting times by ethnicity
colon_cohort %>%
  group_by(ethnicity_group_broad) %>%
  summarise(
    n             = n(),
    median_dx_tx  = median(wt_dx_to_tx,  na.rm = TRUE),
    pct_over_62   = round(100 * mean(wt_dx_to_tx > 62,  na.rm = TRUE), 1),
    pct_over_104  = round(100 * mean(wt_dx_to_tx > 104, na.rm = TRUE), 1)
  ) %>%
  arrange(desc(n)) %>%
  print()

# -----------------------------------------------------------------------------
# 4. Route to diagnosis
# -----------------------------------------------------------------------------

cat("\n=== ROUTE TO DIAGNOSIS ===\n")
colon_cohort %>%
  count(route_combined) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n)) %>%
  print()

# Route by year - detect COVID-era shifts
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  count(year_diag, route_combined) %>%
  group_by(year_diag) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(year_diag, route_combined) %>%
  print(n = Inf)

# Waiting times by route
colon_cohort %>%
  group_by(route_combined) %>%
  summarise(
    n             = n(),
    median_dx_dtt = median(wt_dx_to_dtt, na.rm = TRUE),
    median_dtt_tx = median(wt_dtt_to_tx, na.rm = TRUE),
    median_dx_tx  = median(wt_dx_to_tx,  na.rm = TRUE),
    pct_over_62   = round(100 * mean(wt_dx_to_tx > 62,  na.rm = TRUE), 1),
    pct_over_104  = round(100 * mean(wt_dx_to_tx > 104, na.rm = TRUE), 1)
  ) %>%
  arrange(desc(n)) %>%
  print()

# -----------------------------------------------------------------------------
# 5. Socioeconomic deprivation (IMD quintile)
# -----------------------------------------------------------------------------

cat("\n=== IMD QUINTILE ===\n")

# Note: NHSE_reversed_imd_quintile_lsoas - higher = less deprived
colon_cohort %>%
  count(NHSE_reversed_imd_quintile_lsoas) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(NHSE_reversed_imd_quintile_lsoas) %>%
  print()

# Missing IMD
cat("Missing IMD quintile:",
    sum(is.na(colon_cohort$NHSE_reversed_imd_quintile_lsoas)), "\n")

# Waiting times by IMD quintile
colon_cohort %>%
  group_by(NHSE_reversed_imd_quintile_lsoas) %>%
  summarise(
    n             = n(),
    median_dx_dtt = median(wt_dx_to_dtt, na.rm = TRUE),
    median_dtt_tx = median(wt_dtt_to_tx, na.rm = TRUE),
    median_dx_tx  = median(wt_dx_to_tx,  na.rm = TRUE),
    pct_over_62   = round(100 * mean(wt_dx_to_tx > 62,  na.rm = TRUE), 1),
    pct_over_104  = round(100 * mean(wt_dx_to_tx > 104, na.rm = TRUE), 1)
  ) %>%
  arrange(NHSE_reversed_imd_quintile_lsoas) %>%
  print()

# IMD by year - check for compositional shift
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  group_by(year_diag, NHSE_reversed_imd_quintile_lsoas) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(year_diag) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(year_diag, NHSE_reversed_imd_quintile_lsoas) %>%
  print(n = Inf)

colon_cohort %>%
  mutate(
    year_diag = as.integer(format(diagmdy, "%Y")),
    period    = case_when(
      year_diag <= 2019 ~ "Pre-COVID (2015-19)",
      year_diag == 2020 ~ "COVID (2020)",
      year_diag >= 2021 ~ "Post-COVID (2021-22)"
    )
  ) %>%
  group_by(period, NHSE_reversed_imd_quintile_lsoas) %>%
  summarise(
    n             = n(),
    median_dx_dtt = median(wt_dx_to_dtt, na.rm = TRUE),
    median_dtt_tx = median(wt_dtt_to_tx, na.rm = TRUE),
    median_dx_tx  = median(wt_dx_to_tx,  na.rm = TRUE),
    pct_over_62   = round(100 * mean(wt_dx_to_tx > 62,  na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(period, NHSE_reversed_imd_quintile_lsoas) %>%
  print(n = Inf)


# -----------------------------------------------------------------------------
# 6. Screening detection
# -----------------------------------------------------------------------------

cat("\n=== SCREENING DETECTION ===\n")
colon_cohort %>%
  count(screendetected) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

# Screening by year - bowel cancer screening programme expansion
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  count(year_diag, screendetected) %>%
  group_by(year_diag) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(year_diag, screendetected) %>%
  print(n = Inf)

# Waiting times for screen-detected vs symptomatic
colon_cohort %>%
  group_by(screendetected) %>%
  summarise(
    n             = n(),
    median_dx_dtt = median(wt_dx_to_dtt, na.rm = TRUE),
    median_dtt_tx = median(wt_dtt_to_tx, na.rm = TRUE),
    median_dx_tx  = median(wt_dx_to_tx,  na.rm = TRUE),
    pct_over_62   = round(100 * mean(wt_dx_to_tx > 62,  na.rm = TRUE), 1),
    pct_over_104  = round(100 * mean(wt_dx_to_tx > 104, na.rm = TRUE), 1)
  ) %>%
  print()

# -----------------------------------------------------------------------------
# 7. Stage distribution
# -----------------------------------------------------------------------------

cat("\n=== STAGE DISTRIBUTION ===\n")
colon_cohort %>%
  count(stage) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

# Stage by year - check for COVID-era stage shift
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  count(year_diag, stage) %>%
  group_by(year_diag) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(year_diag, stage) %>%
  print()

# Waiting times by stage
colon_cohort %>%
  group_by(stage) %>%
  summarise(
    n             = n(),
    median_dx_dtt = median(wt_dx_to_dtt, na.rm = TRUE),
    median_dtt_tx = median(wt_dtt_to_tx, na.rm = TRUE),
    median_dx_tx  = median(wt_dx_to_tx,  na.rm = TRUE),
    pct_over_62   = round(100 * mean(wt_dx_to_tx > 62,  na.rm = TRUE), 1),
    pct_over_104  = round(100 * mean(wt_dx_to_tx > 104, na.rm = TRUE), 1)
  ) %>%
  print()

# =============================================================================
# WAITING TIMES BY YEAR - stratified by key characteristics
# =============================================================================

wt_by_year <- function(data, group_var) {
  data %>%
    mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
    group_by(year_diag, {{ group_var }}) %>%
    summarise(
      n             = n(),
      median_dx_dtt = median(wt_dx_to_dtt, na.rm = TRUE),
      median_dtt_tx = median(wt_dtt_to_tx, na.rm = TRUE),
      median_dx_tx  = median(wt_dx_to_tx,  na.rm = TRUE),
      pct_over_62   = round(100 * mean(wt_dx_to_tx > 62,  na.rm = TRUE), 1),
      pct_over_104  = round(100 * mean(wt_dx_to_tx > 104, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    arrange(year_diag, {{ group_var }})
}

# Overall
cat("=== OVERALL ===\n")
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  group_by(year_diag) %>%
  summarise(
    n             = n(),
    median_dx_dtt = median(wt_dx_to_dtt, na.rm = TRUE),
    median_dtt_tx = median(wt_dtt_to_tx, na.rm = TRUE),
    median_dx_tx  = median(wt_dx_to_tx,  na.rm = TRUE),
    pct_over_62   = round(100 * mean(wt_dx_to_tx > 62,  na.rm = TRUE), 1),
    pct_over_104  = round(100 * mean(wt_dx_to_tx > 104, na.rm = TRUE), 1)
  ) %>%
  print()

# By sex
cat("\n=== BY SEX ===\n")
wt_by_year(colon_cohort, sex) %>% print()

# By age group
cat("\n=== BY AGE GROUP ===\n")
colon_cohort %>%
  mutate(age_group = cut(agediag,
                         breaks = c(0, 49, 59, 69, 79, Inf),
                         labels = c("<50", "50-59", "60-69", "70-79", "80+"),
                         right  = TRUE)) %>%
  wt_by_year(age_group) %>%
  print(n = Inf)

# By stage
cat("\n=== BY STAGE ===\n")
wt_by_year(colon_cohort, stage) %>% print()

# By route
cat("\n=== BY ROUTE ===\n")
wt_by_year(colon_cohort, route_combined) %>% print(n = Inf)

# By IMD quintile
cat("\n=== BY IMD QUINTILE ===\n")
wt_by_year(colon_cohort, NHSE_reversed_imd_quintile_lsoas) %>% print(n = Inf)

# By screening detection
cat("\n=== BY SCREENING DETECTION ===\n")
wt_by_year(colon_cohort, screendetected) %>% print()

# By procedure type
cat("\n=== BY PROCEDURE TYPE ===\n")
wt_by_year(colon_cohort, colon_proc_type) %>% print(n = Inf)


