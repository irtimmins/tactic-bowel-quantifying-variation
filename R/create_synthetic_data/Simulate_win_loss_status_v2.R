library(tidyverse)
library(haven)

set.seed(100)

# ── Helper functions ──────────────────────────────────────────────────────────
poisson_pval <- function(n_arrive, n_leave) {
  total <- n_arrive + n_leave
  if (total == 0) return(NA_real_)
  if (n_arrive >= n_leave) {
    pbinom(n_arrive - 1, size = total, prob = 0.5, lower.tail = FALSE)
  } else {
    pbinom(n_arrive,     size = total, prob = 0.5, lower.tail = TRUE)
  }
}

tabulate_wt <- function(df, group_var) {
  df %>%
    filter(!is.na(.data[[group_var]])) %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n = n(),
      across(
        all_of(outcomes),
        list(
          median = ~ median(.x, na.rm = TRUE),
          q25    = ~ quantile(.x, 0.25, na.rm = TRUE),
          q75    = ~ quantile(.x, 0.75, na.rm = TRUE)
        ),
        .names = "{.col}__{.fn}"
      ),
      .groups = "drop"
    ) %>%
    rename(group = 1) %>%
    pivot_longer(-c(group, n), names_to = c("outcome", "stat"), names_sep = "__") %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(result = sprintf("%.0f (%.0f–%.0f)", median, q25, q75)) %>%
    select(group, n, outcome, result) %>%
    pivot_wider(names_from = outcome, values_from = result)
}

kw_pvals <- function(df, group_var) {
  map_dfr(outcomes, function(oc) {
    kt <- kruskal.test(
      reformulate(group_var, response = oc),
      data = df %>% filter(!is.na(.data[[group_var]]))
    )
    tibble(outcome = oc, p.value = round(kt$p.value, 3))
  })
}

# ── Parameters ────────────────────────────────────────────────────────────────
n_sites    <- 100
n_lsoas    <- 32843
n_patients <- 25000

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: SIMULATE SITES AND LSOAS
# ══════════════════════════════════════════════════════════════════════════════

sites <- tibble(
  site_code = paste0("SITE", str_pad(1:n_sites, 3, pad = "0"))
)

lsoa_codes <- paste0("E0100", str_pad(1:n_lsoas, 5, pad = "0"))

lsoas <- tibble(
  lsoa11_code  = lsoa_codes,
  nearest_site = sample(sites$site_code, n_lsoas, replace = TRUE)
)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: SIMULATE DISTANCE MATRIX
# Each LSOA has 5-20 accessible sites, nearest site always included
# Tiny unique increment added per row to break floating point ties
# ensuring R and Stata always identify the same nearest site
# ══════════════════════════════════════════════════════════════════════════════

sim_distance_matrix <- lsoas %>%
  rowwise() %>%
  mutate(
    n_accessible     = sample(5:20, 1),
    other_sites      = list(sample(
      sites$site_code[sites$site_code != nearest_site],
      n_accessible - 1
    )),
    accessible_sites = list(c(nearest_site, other_sites))
  ) %>%
  ungroup() %>%
  select(lsoa11_code, nearest_site, accessible_sites) %>%
  unnest(accessible_sites) %>%
  rename(sitecode = accessible_sites) %>%
  mutate(
    is_nearest       = sitecode == nearest_site,
    total_drive_time = if_else(
      is_nearest,
      runif(n(), min = 1,  max = 15),
      runif(n(), min = 10, max = 90)
    )
  ) %>%
  group_by(lsoa11_code) %>%
  mutate(
    min_non_nearest  = min(total_drive_time[!is_nearest]),
    total_drive_time = if_else(
      is_nearest,
      pmin(total_drive_time, min_non_nearest - runif(1, 1, 5)),
      total_drive_time
    ),
    # Break ties: add unique row increment so no two sites share same drive time
    row_rank         = row_number(),
    total_drive_time = round(total_drive_time, 1) + row_rank * 0.0001
  ) %>%
  ungroup() %>%
  select(lsoa11_code, sitecode, total_drive_time)

# Verify no ties
stopifnot(
  sim_distance_matrix %>%
    group_by(lsoa11_code) %>%
    summarise(n_ties = sum(duplicated(total_drive_time)), .groups = "drop") %>%
    pull(n_ties) %>%
    sum() == 0
)
cat("Distance matrix: no ties confirmed\n")

# Pre-compute accessible sites per LSOA for patient simulation
lsoa_accessible <- sim_distance_matrix %>%
  group_by(lsoa11_code) %>%
  summarise(accessible = list(sitecode), .groups = "drop")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: SIMULATE PATIENTS
# 75% go to nearest site, 25% bypass to a random accessible site
# diag_hosp is always guaranteed to be in the distance matrix
# ══════════════════════════════════════════════════════════════════════════════

patients <- tibble(
  pseudo_patientid = paste0("PAT", str_pad(1:n_patients, 6, pad = "0")),
  lsoa11_code      = sample(lsoa_codes, n_patients, replace = TRUE),
  age              = round(rnorm(n_patients, mean = 67, sd = 11)),
  sex              = sample(c("Male","Female"), n_patients,
                            replace = TRUE, prob = c(0.65, 0.35)),
  stage            = sample(1:3, n_patients, replace = TRUE,
                            prob = c(0.25, 0.40, 0.35)),
  imd_quintile     = sample(1:5, n_patients, replace = TRUE),
  ydiag            = sample(2016:2018, n_patients, replace = TRUE)
) %>%
  left_join(lsoas,           by = "lsoa11_code") %>%
  left_join(lsoa_accessible, by = "lsoa11_code") %>%
  rowwise() %>%
  mutate(
    bypass    = rbinom(1, 1, prob = 0.25),
    diag_hosp = if_else(
      bypass == 0,
      nearest_site,
      sample(setdiff(accessible, nearest_site), 1)
    )
  ) %>%
  ungroup() %>%
  select(-nearest_site, -accessible, -bypass)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: COMPUTE WIN/LOSS FROM COHORT + DISTANCE MATRIX
# ══════════════════════════════════════════════════════════════════════════════

# Nearest site per patient from distance matrix
nearest_from_dm <- sim_distance_matrix %>%
  group_by(lsoa11_code) %>%
  slice_min(total_drive_time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(lsoa11_code, nearest_site_dm = sitecode)

patients <- patients %>%
  left_join(nearest_from_dm, by = "lsoa11_code")

# Sanity checks
stopifnot(sum(is.na(patients$nearest_site_dm)) == 0)

total_balance <- patients %>%
  summarise(
    total_leavers  = sum(diag_hosp != nearest_site_dm),
    total_arrivers = sum(diag_hosp != nearest_site_dm)
  )
cat("Leavers == Arrivers:", total_balance$total_leavers == total_balance$total_arrivers, "\n")

# Aggregate to site level
site_level <- patients %>%
  group_by(nearest_site_dm) %>%
  summarise(n_leavers = sum(diag_hosp != nearest_site_dm), .groups = "drop") %>%
  rename(site_code = nearest_site_dm) %>%
  full_join(
    patients %>%
      group_by(diag_hosp) %>%
      summarise(n_arrivers = sum(diag_hosp != nearest_site_dm), .groups = "drop") %>%
      rename(site_code = diag_hosp),
    by = "site_code"
  ) %>%
  mutate(
    n_leavers  = replace_na(n_leavers,  0),
    n_arrivers = replace_na(n_arrivers, 0),
    n_net_gain = n_arrivers - n_leavers
  ) %>%
  rowwise() %>%
  mutate(p_value = poisson_pval(n_arrivers, n_leavers)) %>%
  ungroup() %>%
  mutate(
    competitor_status = case_when(
      n_net_gain >  0 & p_value <= 0.05 ~ "Winner",
      n_net_gain <  0 & p_value <= 0.05 ~ "Loser",
      TRUE                               ~ "Insignificant diff."
    ),
    competitor_status = factor(
      competitor_status,
      levels = c("Winner", "Loser", "Insignificant diff.")
    )
  )

cat("\nR site-level classification:\n")
print(table(site_level$competitor_status))
site_level %>%
  summarise(
    total_leavers  = sum(n_leavers),
    total_arrivers = sum(n_arrivers),
    n_winner       = sum(competitor_status == "Winner"),
    n_loser        = sum(competitor_status == "Loser"),
    n_insig        = sum(competitor_status == "Insignificant diff.")
  ) %>%
  print()

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: SIMULATE WAITING TIMES BASED ON COMPETITOR STATUS
# Winners have shorter waiting times, Losers longer
# ══════════════════════════════════════════════════════════════════════════════

patients <- patients %>%
  left_join(site_level %>% select(site_code, competitor_status),
            by = c("diag_hosp" = "site_code")) %>%
  mutate(
    wt_mean_offset = case_when(
      competitor_status == "Winner"              ~  -4,
      competitor_status == "Loser"               ~   4,
      competitor_status == "Insignificant diff." ~   0,
      TRUE                                       ~   0
    ),
    wt_dx_to_dtt = pmax(0, round(rnorm(n(), mean = 22 + wt_mean_offset, sd = 10))),
    wt_dtt_to_tx = pmax(0, round(rnorm(n(), mean = 18 + wt_mean_offset, sd = 8))),
    wt_dx_to_tx  = wt_dx_to_dtt + wt_dtt_to_tx
  ) %>%
  filter(
    wt_dx_to_tx  > 0, wt_dx_to_tx <= 180,
    wt_dx_to_dtt >= 0,
    wt_dtt_to_tx >= 0
  ) %>%
  select(-wt_mean_offset)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: TABULATE
# ══════════════════════════════════════════════════════════════════════════════

outcomes <- c("wt_dx_to_tx", "wt_dx_to_dtt", "wt_dtt_to_tx")

cat("\nWaiting times by competitor status:\n")
print(tabulate_wt(patients, "competitor_status"))
print(kw_pvals(patients,    "competitor_status"))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7: EXPORT FOR STATA
# Export cohort WITHOUT competitor_status - Stata derives it independently
# ══════════════════════════════════════════════════════════════════════════════

sim_cohort <- patients %>%
  select(pseudo_patientid, lsoa11_code, diag_hosp,
         age, sex, stage, imd_quintile, ydiag,
         wt_dx_to_dtt, wt_dtt_to_tx, wt_dx_to_tx)

write_dta(sim_cohort,
          "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/provider_level/simulated_colon_cohort.dta")
write_dta(sim_distance_matrix,
          "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/provider_level/simulated_distance_matrix.dta")

cat("\nFiles exported. Run Stata script to verify identical classification.\n")


site_level %>%
  summarise(
    total_leavers  = sum(n_leavers),
    total_arrivers = sum(n_arrivers),
    n_winner       = sum(competitor_status == "Winner"),
    n_loser        = sum(competitor_status == "Loser"),
    n_insig        = sum(competitor_status == "Insignificant diff.")
  )
