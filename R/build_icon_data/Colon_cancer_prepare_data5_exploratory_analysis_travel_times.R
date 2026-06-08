
library(tidyverse)
library(haven)

distance_matrix_v1 <- read_dta(
   "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/provider_level/bowel_pairwise_distance_matrix.dta"
) %>%
   rename(lsoa11_code = lsoa11cd)

# #head(distance_matrix)

library(readr)
distance_matrix_v2 <- read_csv("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/provider_level/SORT drive time matrix_all hospitals_07.08.2024.csv") %>%
  rename(lsoa11_code = lsoa11cd ) %>%
  rename(sitecode =hospital_code)
distance_matrix_v2
names(distance_matrix_v1)
names(distance_matrix_v2)
distance_matrix_v2 %>%
  mutate(n_decimals = nchar(sub(".*\\.", "", as.character(total_drive_time)))) %>%
  count(n_decimals)

head(distance_matrix_v1)
head(distance_matrix_v2)

list1 <- unique(distance_matrix_v1$sitecode)
list2 <- unique(distance_matrix_v2$sitecode)
sum(is.element(list1, list2))
sum(is.element(list2, list1))


# ── Coverage ───────────────────────────────────────────────────────────────────
cat("V1 sites:", length(list1), "\n")  
cat("V2 sites:", length(list2), "\n")
cat("Shared sites:", sum(is.element(list1, list2)), "\n")
cat("V1 only:", paste(setdiff(list1, list2), collapse = ", "), "\n")
cat("V2 only:", paste(setdiff(list2, list1), collapse = ", "), "\n")

# ── Drive time concordance on shared LSOA x site pairs ────────────────────────
joined <- distance_matrix_v1 %>%
  select(lsoa11_code, sitecode, dt_v1 = total_drive_time) %>%
  inner_join(
    distance_matrix_v2 %>% select(lsoa11_code, sitecode, dt_v2 = total_drive_time),
    by = c("lsoa11_code", "sitecode")
  ) %>%
  mutate(diff = dt_v2 - dt_v1)

cat("Matched pairs:", nrow(joined), "\n")
cat("V1 only pairs:", nrow(distance_matrix_v1) - nrow(joined), "\n")
cat("V2 only pairs:", nrow(distance_matrix_v2) - nrow(joined), "\n")

summary(joined$diff)
cor(joined$dt_v1, joined$dt_v2, method = "spearman")

# ── Does nearest hospital agree? ──────────────────────────────────────────────
nearest_v1 <- distance_matrix_v1 %>%
  group_by(lsoa11_code) %>%
  slice_min(total_drive_time, n = 1, with_ties = FALSE) %>%
  select(lsoa11_code, nearest_v1 = sitecode)

nearest_v2 <- distance_matrix_v2 %>%
  group_by(lsoa11_code) %>%
  slice_min(total_drive_time, n = 1, with_ties = FALSE) %>%
  select(lsoa11_code, nearest_v2 = sitecode)

nearest_v1 %>%
  inner_join(nearest_v2, by = "lsoa11_code") %>%
  ungroup() %>%
  summarise(
    n         = n(),
    pct_agree = mean(nearest_v1 == nearest_v2) * 100
  )


# Where do they disagree - how different are the drive times?
nearest_v1 %>%
  inner_join(nearest_v2, by = "lsoa11_code") %>%
  ungroup() %>%
  filter(nearest_v1 != nearest_v2) %>%
  inner_join(distance_matrix_v1 %>% select(lsoa11_code, sitecode, dt_v1 = total_drive_time),
             by = c("lsoa11_code", "nearest_v1" = "sitecode")) %>%
  inner_join(distance_matrix_v2 %>% select(lsoa11_code, sitecode, dt_v2 = total_drive_time),
             by = c("lsoa11_code", "nearest_v2" = "sitecode")) %>%
  mutate(dt_diff = abs(dt_v1 - dt_v2)) %>%
  summarise(
    n_disagree    = n(),
    median_diff   = median(dt_diff),
    mean_diff     = mean(dt_diff),
    pct_within_5  = mean(dt_diff <= 5) * 100,
    pct_within_10 = mean(dt_diff <= 10) * 100
  )


# ── Load data ──────────────────────────────────────────────────────────────────
colon_cohort <- readRDS(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_cci_2015_2022.rds"
) %>%
  filter(ydiag %in% 2020:2022) 
# 

distance_matrix <- read_dta(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/provider_level/bowel_pairwise_distance_matrix.dta"
) %>%
  rename(lsoa11_code = lsoa11cd)


# ── Merge cohort with distance matrix ─────────────────────────────────────────
# Each patient (identified by lsoa11_code) is matched to all hospitals
# with their drive times — this creates the long-format patient x site data

# Hospital inclusion
hosp_min <- colon_cohort  %>%
  count(diag_hosp, ydiag) %>%
  group_by(diag_hosp) %>%
  summarise(
    n_years    = n_distinct(ydiag),
    min_annual = min(n),
    .groups    = "drop"
  ) %>%
  filter(
    min_annual >= 5
  ) %>%
  pull(diag_hosp)

cat("Hospitals meeting inclusion criteria:", length(hosp_min), "\n")

colon_cohort <- colon_cohort %>%
  filter(diag_hosp %in% hosp_min) 

hosp_list <- unique(distance_matrix$sitecode)
hosp_list_cohort <- unique(colon_cohort$diag_hosp)
sum(is.element(hosp_list_cohort, hosp_list))

df <- colon_cohort %>%
  select(pseudo_patientid, lsoa11_code, diag_hosp) %>%   # SITETRET = treating hospital
  inner_join(distance_matrix, by = "lsoa11_code", relationship = "many-to-many") %>%
  rename(
    site_code      = sitecode,
    drive_time     = total_drive_time,
    site_diag     = diag_hosp        # 1 if this site_code is where diagnosis occurred
  ) %>%
  mutate(site_resec = as.integer(site_code == site_diag)) %>%
  filter(site_code %in% hosp_list_cohort)

#head(df)

# ── Identify nearest site per patient ─────────────────────────────────────────
df <- df %>%
  group_by(pseudo_patientid) %>%
  mutate(
    drive_time_min = min(drive_time),
    site_nearest   = as.integer(drive_time == drive_time_min)
  ) %>%
  ungroup()

# ── Core patients count per site ──────────────────────────────────────────────
site_core <- df %>%
  filter(site_nearest == 1) %>%
  group_by(site_code) %>%
  summarise(npat_core = n(), .groups = "drop")

# ── Keep only relevant rows (leavers + arrivers + core treated) ───────────────
# Drop rows where patient neither had resection at this site nor is nearest site
df <- df %>%
  filter(!(site_resec == 0 & site_nearest == 0))

# ── Leavers and arrivers ───────────────────────────────────────────────────────
df <- df %>%
  mutate(
    leave  = as.integer(site_resec == 0 & site_nearest == 1),
    arrive = as.integer(site_resec == 1 & site_nearest == 0)
  )

# ── Aggregate to site level ────────────────────────────────────────────────────
site_level <- df %>%
  group_by(site_code) %>%
  summarise(
    n_leave  = sum(leave),
    n_arrive = sum(arrive),
    .groups  = "drop"
  ) %>%
  left_join(site_core, by = "site_code") %>%
  mutate(n_net = n_arrive - n_leave)

site_level


# Do we have any negative n_net at all?
summary(site_level$n_net)

# How many sites have n_net < 0?
sum(site_level$n_net < 0)

# ── Test for significant difference: arrivers vs leavers (1-sided Poisson) ────
# Equivalent to Stata's -iri- conditional test (Krishnamoorthy & Thomson 2004)
poisson_pval <- function(n_arrive, n_leave) {
  total <- n_arrive + n_leave
  if (total == 0) return(NA_real_)
  
  # One-sided test in the appropriate direction
  if (n_arrive >= n_leave) {
    # Test if arrivers significantly greater than leavers
    pbinom(n_arrive - 1, size = total, prob = 0.5, lower.tail = FALSE)
  } else {
    # Test if leavers significantly greater than arrivers
    pbinom(n_arrive, size = total, prob = 0.5, lower.tail = TRUE)
  }
}

# Recalculate
site_level <- site_level %>%
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

table(site_level$win_lose)

site_level <- site_level %>%
  rowwise() %>%
  mutate(p_value = poisson_pval(n_arrive, n_leave)) %>%
  ungroup()

# ── Classify winner / loser / insignificant ────────────────────────────────────
site_level <- site_level %>%
  mutate(
    win_lose = case_when(
      n_net >  0 & p_value <= 0.05 ~ "Winner",
      n_net <  0 & p_value <= 0.05 ~ "Loser",
      TRUE                          ~ "Insignificant diff."
    ),
    win_lose = factor(win_lose, levels = c("Winner", "Loser", "Insignificant diff."))
  )

# ── Check ──────────────────────────────────────────────────────────────────────
cat("Number of sites:", nrow(site_level), "\n")
cat("Sites with n_net == 0:", sum(site_level$n_net == 0), "\n")
print(table(site_level$win_lose, useNA = "always"))


# ── Save ───────────────────────────────────────────────────────────────────────
#saveRDS(site_level, "SCI_colon.rds")

names(df)
site_level
cohort_df <- colon_cohort %>%
  filter(
    !is.na(wt_dx_to_tx), !is.na(wt_dx_to_dtt), !is.na(wt_dtt_to_tx),
    wt_dx_to_tx  > 0, wt_dx_to_tx <= 180,
    wt_dx_to_dtt >= 0,
    wt_dtt_to_tx >= 0
  ) %>%
  mutate(
    stage          = factor(stage, levels = 1:3,
                            labels = c("Stage I","Stage II","Stage III")),
    route_combined = factor(route_combined),
    sex            = factor(sex, labels = c("Male","Female")),
    imd_q          = factor(NHSE_reversed_imd_quintile_lsoas,
                            levels = c("1 - most deprived","2","3","4","5 - least deprived"),
                            labels = c("Q1 Most\ndeprived","Q2","Q3","Q4","Q5 Least\ndeprived")),
    ydiag          = as.integer(ydiag),
    age_group      = cut(agediag, breaks = c(0,49,59,69,79,Inf),
                         labels = c("<50","50-59","60-69","70-79","80+"),
                         right  = TRUE)
  )
cohort_df$wt_dx_to_tx

##############################################################

#install.packages("gtsummary")
#install.packages("gt")
library(gtsummary)
library(gt)

# ── Join site-level factors to patient cohort ──────────────────────────────────
cohort_df <- cohort_df %>%
  left_join(
    site_level %>% select(site_code, win_lose, n_net),
    by = c("SITETRET" = "site_code")
  )

# ── Outcomes to model ──────────────────────────────────────────────────────────
outcomes <- c("wt_dx_to_tx", "wt_dx_to_dtt", "wt_dtt_to_tx")

# ── Covariates ─────────────────────────────────────────────────────────────────
covars <- c(
  "age_group", "sex", "imd_q", "stage",
  "route_combined", "ydiag", "cci_group"
)


# ── Helper: median (IQR) by group ──────────────────────────────────────────────
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

# ── Call with strings ──────────────────────────────────────────────────────────
tab_winlose <- tabulate_wt(cohort_df, "win_lose")
pvals_winlose <- kw_pvals(cohort_df, "win_lose")

# ── Tabulate by winner/loser ───────────────────────────────────────────────────
tab_winlose <- tabulate_wt(cohort_df, "win_lose")
print(tab_winlose)

# ── Kruskal-Wallis p-values ────────────────────────────────────────────────────
kw_pvals <- function(df, group_var) {
  map_dfr(outcomes, function(oc) {
    kt <- kruskal.test(
      x = df[[oc]],
      g = df[[group_var]]
    )
    tibble(
      outcome = oc,
      p.value = kt$p.value #, 
      # p.value = case_when(
      #   kt$p.value < 0.001 ~ "<0.001",
      #   kt$p.value < 0.01  ~ as.character(round(kt$p.value, 3)),
      #   TRUE               ~ as.character(round(kt$p.value, 2))
      #)
    )
  })
}

pvals_winlose <- kw_pvals(cohort_df, "win_lose")
pvals_winlose


# ── Pretty gt tables ───────────────────────────────────────────────────────────
format_gt <- function(tab, pvals, title) {
  
  pval_row <- pvals %>%
    pivot_wider(names_from = outcome, values_from = p.value) %>%
    mutate(
      group = "p-value", 
      n = NA_integer_,
      across(any_of(outcomes), ~ as.character(.x))
    )
  
  tab %>%
    mutate(across(any_of(outcomes), as.character)) %>%
    bind_rows(pval_row) %>%
    gt() %>%
    tab_header(title = title) %>%
    tab_spanner(
      label   = "Diagnosis to treatment (days)",
      columns = contains("wt_dx_to_tx")
    ) %>%
    tab_spanner(
      label   = "Diagnosis to decision to treat (days)",
      columns = contains("wt_dx_to_dtt")
    ) %>%
    tab_spanner(
      label   = "Decision to treat to treatment (days)",
      columns = contains("wt_dtt_to_tx")
    ) %>%
    cols_label(
      group = "",
      n     = "N"
    ) %>%
    tab_footnote("Values are median (IQR). P-value from Kruskal-Wallis test.")
}

gt_winlose <- format_gt(tab_winlose, pvals_winlose, "Waiting times by competitive status")
gt_winlose
