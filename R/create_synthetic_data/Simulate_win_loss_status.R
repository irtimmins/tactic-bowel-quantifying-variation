library(tidyverse)

set.seed(100)

softmax <- function(x) exp(x) / sum(exp(x))

poisson_pval <- function(n_arrive, n_leave) {
  total <- n_arrive + n_leave
  if (total == 0) return(NA_real_)
  if (n_arrive >= n_leave) {
    pbinom(n_arrive - 1, size = total, prob = 0.5, lower.tail = FALSE)
  } else {
    pbinom(n_arrive, size = total, prob = 0.5, lower.tail = TRUE)
  }
}

# ── Helper: median (IQR) by group ──────────────────────────────────────────────
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
    pivot_longer(
      cols      = -c(group, n),
      names_to  = c("outcome", "stat"),
      names_sep = "__"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(result = sprintf("%.0f (%.0f–%.0f)", median, q25, q75)) %>%
    select(group, n, outcome, result) %>%
    pivot_wider(names_from = outcome, values_from = result)
}

# ── Kruskal-Wallis p-values ────────────────────────────────────────────────────
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

# ── Simulate sites ────────────────────────────────────────────────────────────
sites <- tibble(
  site_code      = paste0("SITE", str_pad(1:n_sites, 3, pad = "0")),
  attractiveness = rnorm(n_sites, mean = 0, sd = 1)
)

# ── Simulate LSOAs ────────────────────────────────────────────────────────────

# Generate LSOAs once
lsoa_codes <- paste0("E0100", str_pad(1:n_lsoas, 5, pad = "0"))
# gives "E010000001" to "E010032843" - consistent across both files

lsoas <- tibble(
  lsoa11_code  = lsoa_codes,
  nearest_site = sample(sites$site_code, n_lsoas, replace = TRUE)
)
# ── Simulate distance matrix ──────────────────────────────────────────────────
sim_distance_matrix <- lsoas %>%
  rowwise() %>%
  mutate(
    n_sites_accessible = sample(5:20, 1),
    # Always include nearest_site, then sample the rest
    other_sites        = list(sample(
      sites$site_code[sites$site_code != nearest_site], 
      n_sites_accessible - 1
    )),
    accessible_sites   = list(c(nearest_site, other_sites))
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
    )
  ) %>%
  ungroup() %>%
  mutate(total_drive_time = pmax(0.5, round(total_drive_time, 1))) %>%
  select(lsoa11_code, sitecode, total_drive_time)

# ── Simulate patients ─────────────────────────────────────────────────────────
patients <- tibble(
  pseudo_patientid = paste0("PAT", str_pad(1:n_patients, 6, pad = "0")),
  lsoa11_code      = sample(lsoa_codes, n_patients, replace = TRUE),
  age              = round(rnorm(n_patients, mean = 67, sd = 11)),
  sex              = sample(c("Male","Female"), n_patients, replace = TRUE, prob = c(0.65, 0.35)),
  stage            = sample(1:3, n_patients, replace = TRUE, prob = c(0.25, 0.40, 0.35)),
  imd_quintile     = sample(1:5, n_patients, replace = TRUE),
  ydiag            = sample(2016:2018, n_patients, replace = TRUE)
) %>%
  left_join(lsoas, by = "lsoa11_code")

# ── Assign diagnosing hospital based on attractiveness ────────────────────────
# Pre-compute which sites are accessible from each LSOA
lsoa_accessible_sites <- sim_distance_matrix %>%
  group_by(lsoa11_code) %>%
  summarise(accessible = list(sitecode), .groups = "drop")

patients <- patients %>%
  left_join(sites %>% rename(nearest_site    = site_code,
                             attract_nearest = attractiveness),
            by = "nearest_site") %>%
  left_join(lsoa_accessible_sites, by = "lsoa11_code") %>%
  rowwise() %>%
  mutate(
    bypass    = rbinom(1, 1, prob = plogis(-attract_nearest * 0.5 + 0.3)),
    # When bypassing, only sample from sites in this LSOA's distance matrix
    diag_hosp = if_else(
      bypass == 0,
      nearest_site,
      {
        candidates <- setdiff(accessible, nearest_site)
        candidate_attract <- sites$attractiveness[sites$site_code %in% candidates]
        sample(candidates, 1, prob = softmax(candidate_attract))
      }
    )
  ) %>%
  ungroup() %>%
  select(-accessible) %>%
  left_join(sites %>% rename(diag_hosp    = site_code,
                             attract_diag = attractiveness),
            by = "diag_hosp")

# ── Simulate waiting times ────────────────────────────────────────────────────
patients <- patients %>%
  mutate(
    wt_dx_to_dtt = pmax(0, round(rnorm(n(), mean = 22 - attract_diag * 2,   sd = 10))),
    wt_dtt_to_tx = pmax(0, round(rnorm(n(), mean = 18 - attract_diag * 1.5, sd = 8))),
    wt_dx_to_tx  = wt_dx_to_dtt + wt_dtt_to_tx
  ) %>%
  filter(
    wt_dx_to_tx  > 0, wt_dx_to_tx <= 180,
    wt_dx_to_dtt >= 0,
    wt_dtt_to_tx >= 0
  )

# ── Derive nearest hospital from distance matrix ──────────────────────────────
nearest_from_dm <- sim_distance_matrix %>%
  group_by(lsoa11_code) %>%
  slice_min(total_drive_time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(lsoa11_code, nearest_site_dm = sitecode)

# ── Join nearest site (from distance matrix) back to patients ─────────────────
patients <- patients %>%
  left_join(nearest_from_dm, by = "lsoa11_code")

# ── Derive winner/loser from distance matrix nearest site ─────────────────────
sim_site_level <- patients %>%
  # leavers: nearest site per DM but treated elsewhere
  group_by(nearest_site_dm) %>%
  summarise(n_leave = sum(diag_hosp != nearest_site_dm), .groups = "drop") %>%
  rename(site_code = nearest_site_dm) %>%
  full_join(
    # arrivers: treated here but nearest site per DM was elsewhere
    patients %>%
      group_by(diag_hosp) %>%
      summarise(n_arrive = sum(diag_hosp != nearest_site_dm), .groups = "drop") %>%
      rename(site_code = diag_hosp),
    by = "site_code"
  ) %>%
  mutate(
    n_leave  = replace_na(n_leave, 0),
    n_arrive = replace_na(n_arrive, 0),
    n_net    = n_arrive - n_leave
  ) %>%
  rowwise() %>%
  mutate(p_value = poisson_pval(n_arrive, n_leave)) %>%
  ungroup() %>%
  mutate(
    win_lose = case_when(
      n_net >  0 & p_value <= 0.05 ~ "Winner",
      n_net <  0 & p_value <= 0.05 ~ "Loser",
      TRUE                          ~ "Insignificant diff."
    ),
    win_lose = factor(win_lose, levels = c("Winner", "Loser", "Insignificant diff."))
  )

print(table(sim_site_level$win_lose))

# ── Final cohort ──────────────────────────────────────────────────────────────
sim_cohort <- patients %>%
  select(pseudo_patientid, lsoa11_code, diag_hosp,
         age, sex, stage, imd_quintile, ydiag,
         wt_dx_to_dtt, wt_dtt_to_tx, wt_dx_to_tx) %>%
  left_join(
    sim_site_level %>% select(site_code, win_lose),
    by = c("diag_hosp" = "site_code")
  )

# ── Tabulate ──────────────────────────────────────────────────────────────────
outcomes <- c("wt_dx_to_tx", "wt_dx_to_dtt", "wt_dtt_to_tx")

tab_winlose   <- tabulate_wt(sim_cohort, "win_lose")
pvals_winlose <- kw_pvals(sim_cohort, "win_lose")

print(tab_winlose)
print(pvals_winlose)

library(haven)

write_dta(sim_cohort,          "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/provider_level/simulated_colon_cohort.dta")
write_dta(sim_distance_matrix, "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/provider_level/simulated_distance_matrix.dta")



# # ── Save ──────────────────────────────────────────────────────────────────────
# saveRDS(sim_cohort,          "sim_cohort.rds")
# saveRDS(sim_site_level,      "sim_site_level.rds")
# saveRDS(sim_distance_matrix, "sim_distance_matrix.rds")