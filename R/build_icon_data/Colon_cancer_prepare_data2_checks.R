colon_cohort <- readRDS(

  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_2015_2022.rds"
)


# -----------------------------------------------------------------------------
# 1. Row counts and linkage at each pipeline stage
# -----------------------------------------------------------------------------

cat("=== COHORT ATTRITION ===\n")
cat("NCRAS colon 2015-2022 (all stages/routes): ",
    nrow(ncras_colon), "\n")
cat("NCRAS after stage 1-3 + elective route:    ",
    n_distinct(ncras_index$pseudo_patientid), "\n")
cat("NCRAS + HES (matched elective resection):  ",
    n_distinct(ncras_hes$pseudo_patientid), "\n")
cat("Final cohort (NCRAS + HES + CWT):          ",
    n_distinct(colon_cohort$pseudo_patientid), "\n")

# -----------------------------------------------------------------------------
# 2. Linkage rates by year
# -----------------------------------------------------------------------------

cat("\n=== LINKAGE RATES BY YEAR ===\n")
ncras_index %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  count(year_diag, name = "n_ncras") %>%
  left_join(
    ncras_hes %>%
      mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
      count(year_diag, name = "n_hes"),
    by = "year_diag"
  ) %>%
  left_join(
    colon_cohort %>%
      mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
      count(year_diag, name = "n_cwt"),
    by = "year_diag"
  ) %>%
  mutate(
    pct_hes_linked = round(100 * n_hes / n_ncras, 1),
    pct_cwt_linked = round(100 * n_cwt / n_ncras, 1)
  ) %>%
  print()

# -----------------------------------------------------------------------------
# 3. Pathway sequence validity
# -----------------------------------------------------------------------------

cat("\n=== PATHWAY SEQUENCE VALIDITY ===\n")
colon_cohort %>%
  summarise(across(ends_with("_ok"), ~ round(mean(.x) * 100, 1))) %>%
  print()

# -----------------------------------------------------------------------------
# 4. Waiting time distributions
# -----------------------------------------------------------------------------

cat("\n=== WAITING TIME DISTRIBUTIONS ===\n")
colon_cohort %>%
 # filter(wt_dx_to_dtt >= 0) %>%
  summarise(across(
    c(wt_dx_to_dtt, wt_dtt_to_tx, wt_dx_to_tx),
    list(
      n_nonmiss = ~ sum(!is.na(.x)),
      p25    = ~ quantile(.x, 0.25, na.rm = TRUE),
      median = ~ median(.x, na.rm = TRUE),
      mean   = ~ round(mean(.x, na.rm = TRUE), 1),
      sd = ~ round(sd(.x), 1),
      p75    = ~ quantile(.x, 0.75, na.rm = TRUE),
      p90    = ~ quantile(.x, 0.90, na.rm = TRUE),
      pct_over_62  = ~ round(100 * mean(.x > 62,  na.rm = TRUE), 1),
      pct_over_104 = ~ round(100 * mean(.x > 104, na.rm = TRUE), 1)
    )
  )) %>%
  pivot_longer(everything(),
               names_to  = c("interval", "stat"),
               names_sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  print()

# -----------------------------------------------------------------------------
# 5. Waiting times by year and procedure
# -----------------------------------------------------------------------------

cat("\n=== WAITING TIMES BY YEAR ===\n")
colon_cohort %>%
#  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  group_by(ydiag) %>%
  summarise(
    n             = n(),
    mean_dx_dtt = mean(wt_dx_to_dtt, na.rm = TRUE),
    mean_dtt_tx = mean(wt_dtt_to_tx, na.rm = TRUE),
    median_dtt_tx = median(wt_dtt_to_tx, na.rm = TRUE),
    mean_dx_tx  = mean(wt_dx_to_tx,  na.rm = TRUE),
    pct_over_62   = round(100 * mean(wt_dx_to_tx > 62,  na.rm = TRUE), 1),
    pct_over_104  = round(100 * mean(wt_dx_to_tx > 104, na.rm = TRUE), 1)
  ) %>%
  print()

colon_cohort %>%
  filter(grepl("^R[A-Z0-9]{4}$", diag_hosp)) %>%
  ggplot(aes(x = reorder(diag_hosp, wt_dtt_to_tx, median),  y = wt_dtt_to_tx)) + 
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  ylim(0,100) +
  geom_hline(yintercept = 31)

summary(colon_cohort$wt_dx_to_dtt)
summary(colon_cohort$wt_dtt_to_tx)



cat("\n=== WAITING TIMES BY PROCEDURE ===\n")
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  group_by(colon_proc_type) %>%
  summarise(
    n             = n(),
    mean_dx_dtt = mean(wt_dx_to_dtt, na.rm = TRUE),
    mean_dtt_tx = mean(wt_dtt_to_tx, na.rm = TRUE),
    mean_dx_tx  = mean(wt_dx_to_tx,  na.rm = TRUE),
    pct_over_62   = round(100 * mean(wt_dx_to_tx > 62,  na.rm = TRUE), 1),
    pct_over_104  = round(100 * mean(wt_dx_to_tx > 104, na.rm = TRUE), 1)
  ) %>%
  print()

# -----------------------------------------------------------------------------
# 6. Missing data rates
# -----------------------------------------------------------------------------

cat("\n=== MISSING DATA BY YEAR ===\n")
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  group_by(year_diag) %>%
  summarise(
    n              = n(),
    pct_miss_dtt   = round(100 * mean(is.na(dtt_date)), 1),
    pct_miss_mdt   = round(100 * mean(is.na(mdt_date)), 1),
    pct_miss_stage = round(100 * mean(is.na(stage)),    1),
    pct_miss_grade = round(100 * mean(is.na(grade)),    1)
  ) %>%
  print()

# -----------------------------------------------------------------------------
# 7. Modality by year - confirm coding transition is handled correctly
# -----------------------------------------------------------------------------

cat("\n=== MODALITY BY YEAR ===\n")
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  count(year_diag, modality) %>%
  group_by(year_diag) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(year_diag, modality) %>%
  print(n = Inf)

# -----------------------------------------------------------------------------
# 8. ICD-10 site match between NCRAS and CWT
# -----------------------------------------------------------------------------

cat("\n=== ICD-10 SITE MATCH ===\n")
colon_cohort %>%
  count(site_match) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

colon_cohort %>%
  filter(site_match == 0) %>%
  count(site_icd10, sitestr, sort = TRUE) %>%
  print()

# -----------------------------------------------------------------------------
# 9. Procedure type distribution by year
# -----------------------------------------------------------------------------

cat("\n=== PROCEDURE TYPE BY YEAR ===\n")
colon_cohort %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  count(year_diag, colon_proc_type) %>%
  group_by(year_diag) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(year_diag, colon_proc_type) %>%
  print(n = Inf)

# -----------------------------------------------------------------------------
# 10. Route to diagnosis by year
# -----------------------------------------------------------------------------

cat("\n=== ROUTE TO DIAGNOSIS BY YEAR ===\n")
colon_cohort %>%
 # filter(ydiag %in% c(2020:2022)) %>%
  group_by(route_combined) %>%
  summarise(mean = mean(wt_dx_to_tx),
            n = n()) %>%
 # mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
 # count(route_combined) %>%
 # group_by(year_diag) %>%
  ungroup() %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(route_combined)# %>%
  print(n = Inf)

  colon_cohort %>%
    # filter(ydiag %in% c(2020:2022)) %>%
    group_by(ydiag) %>%
    summarise(mean = mean(wt_dx_to_tx),
              n = n()) %>%
    # mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
    # count(route_combined) %>%
    # group_by(year_diag) %>%
    ungroup() %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
  print(n = Inf)
  
# -----------------------------------------------------------------------------
# 11. One row per patient check
# -----------------------------------------------------------------------------

cat("\n=== DUPLICATE CHECK ===\n")
dupe_count <- colon_cohort %>%
  count(pseudo_patientid) %>%
  filter(n > 1) %>%
  nrow()
cat("Patients with >1 row:", dupe_count,
    if (dupe_count == 0) "--- PASS" else "--- FAIL", "\n")

# -----------------------------------------------------------------------------
# 12. Date range sanity checks
# -----------------------------------------------------------------------------

cat("\n=== DATE RANGES ===\n")
colon_cohort %>%
  summarise(
    diag_min  = min(diagmdy,   na.rm = TRUE),
    diag_max  = max(diagmdy,   na.rm = TRUE),
    tx_min    = min(tx_date,   na.rm = TRUE),
    tx_max    = max(tx_date,   na.rm = TRUE),
    dtt_min   = min(dtt_date,  na.rm = TRUE),
    dtt_max   = max(dtt_date,  na.rm = TRUE)
  ) %>%
  print()


colon_cohort %>% filter(diag_hosp == "RPY01") %>%
  View()
  nrow()


  
  
  recent_data_by_year <- readRDS("D:/Projects/#2045_ICON_TACTIC/Project3_MSc_ses/Prepare_data/raw_means_hospital_colon_surgery_2023_2025.rds")
  recent_data_by_year
  
  View(recent_data_by_year)
  
  
  summary(as.factor(recent_data_by_year$year_start))
  test4 <- recent_data_by_year %>%
    filter(year_start == "July 2024") %>%
    # filter(hosp %in% test3$hosp) %>%
    select(hosp, margin) %>%
    rename(est=  margin) %>%
    mutate(model = 2023)
  summary(test4$est)
  mean(test4$est, na.rm = T)
  hist(test4$est, breaks = 50)
  
  test4
  
  test_all2 <- test3 %>%
    bind_rows(test4) %>%
    pivot_wider(names_from = model, values_from = c(est))#
  #test_all2
  cor(test_all2$`2020`, test_all2$`2023`)
  
  
  test_all2 %>%
    ggplot(aes(x = `2020`, y=`2023`))+
    geom_point()+
    theme_classic()+
    geom_abline(a = 0, b =1 )
  
  
  cor(test_all2$`2022`, test_all2$`2023`)
  #test_all2
  summary(lm(`2023`~`2022`, data = test_all2))
  
  
  
  df
  
  
test <- readRDS("D:/Projects/#2045_ICON_TACTIC/Project3_MSc_ses/Prepare_data/raw_hosp_year_summary_bowel_2023_2025.rds")


test  <- test %>%
  rename(hosp_diag = hosp_diagnosis,
         ydiag = year_diag) %>%
  group_by(hosp_diag, ydiag) %>%
  filter(row_number()==1) %>%
  ungroup()


length(unique(test$hosp_diag))  
View(test)  
  
test %>% filter(hosp_diag == "R0A07")

colon_test <- colon_cohort_qc %>%
  filter(diag_hosp %in% test$hosp_diag)

length(unique(colon_test$diag_hosp))  





cat("\n=== WAITING TIME DISTRIBUTIONS ===\n")

icon_hosp <- colon_test %>%
  mutate(year_diag = as.integer(format(diagmdy, "%Y"))) %>%
  group_by(ydiag, diag_hosp) %>%
  summarise(
    n             = n(),
    mean_dx_dtt = mean(wt_dx_to_dtt, na.rm = TRUE),
  ) %>%
 ungroup() %>%
 mutate(data = "icon") %>%
 filter(diag_hosp %in% rapid_hosp$diag_hosp)

#names(test)

#hosp_means <-  
rapid_hosp <- test %>%
  rename(wt_dx_to_dtt = mean,
         diag_hosp = hosp_diag) %>%
  group_by(ydiag, diag_hosp) %>%
  summarise(
    n             = obs,
    mean_dx_dtt = mean(wt_dx_to_dtt, na.rm = TRUE),
  )  %>%
  ungroup() %>%
  mutate(data = "sarah_rapid") %>%
  filter(diag_hosp %in% icon_hosp$diag_hosp)

compare_data <- icon_hosp %>%
  bind_rows(rapid_hosp)


compare_data %>%
  arrange(diag_hosp, ydiag)

# 
# 
# test %>%
#   rename(wt_dx_to_dtt = mean) %>%
#   group_by(ydiag) %>%
#   summarise(
#     n             = n(),
#     mean_dx_dtt = mean(wt_dx_to_dtt, na.rm = TRUE),
#   ) %>%
#   print()
#   
#   
  
library(tidyverse)

# ── 1. Split by source ──────────────────────────────────────────────────────
icon_data  <- compare_data %>% filter(data == "icon")
sarah_data <- compare_data %>% filter(data == "sarah_rapid")


library(tidyverse)

compare_data %>%
  group_by(ydiag) %>%
  summarise(mean = mean(mean_dx_dtt),
            n= sum(n))


# ── 1. Find the transition year per hospital ────────────────────────────────
transitions <- compare_data %>%
  group_by(diag_hosp) %>%
  summarise(
    last_icon   = max(ydiag[data == "icon"],        na.rm = TRUE),
    first_sarah = min(ydiag[data == "sarah_rapid"], na.rm = TRUE),
    gap_years   = first_sarah - last_icon - 1,
    .groups = "drop"
  )

transitions
# ── 2. Extract the year before and after handover per hospital ──────────────
handover <- compare_data %>%
  left_join(transitions, by = "diag_hosp") %>%
  filter(ydiag == last_icon | ydiag == first_sarah) %>%
  mutate(period = if_else(data == "icon", "pre", "post")) %>%
  select(diag_hosp, period, n, mean_dx_dtt) %>%
  pivot_wider(names_from = period, values_from = c(n, mean_dx_dtt)) %>%
  mutate(
    n_jump        = n_post - n_pre,
    n_pct_jump    = (n_jump / n_pre) * 100,
    mean_dtt_jump = mean_dx_dtt_post - mean_dx_dtt_pre
  ) %>%
  arrange(desc(abs(mean_dtt_jump)))

# ── 3. Flag suspicious hospitals ───────────────────────────────────────────
handover %>%
  filter(abs(n_pct_jump) > 30 | abs(mean_dtt_jump) > 10) %>%
  print(n = 30)

# ── 4. Plot: waiting time jump at handover across all hospitals ────────────
handover %>%
  ggplot(aes(x = mean_dx_dtt_pre, y = mean_dx_dtt_post)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = mean_dtt_jump), size = 2, alpha = 0.8) +
  scale_colour_gradient2(low = "steelblue", mid = "grey80", high = "red",
                         midpoint = 0, name = "wait jump (days)") +
  labs(
    title    = "Waiting time at source handover (icon → sarah_rapid)",
    subtitle = "Points above line = sarah_rapid records longer waits than icon",
    x = "icon: mean wait (last year)",
    y = "sarah_rapid: mean wait (first year)"
  ) +
  theme_minimal()+
  coord_fixed(xlim = c(0, 60), ylim = c(0, 60))

# ── 5. Overall trend plot for a single hospital (swap out code as needed) ──
compare_data %>%
  filter(diag_hosp == "RHU03") %>%
  ggplot(aes(x = ydiag, y = mean_dx_dtt, colour = data, group = 1)) +
  geom_vline(data = transitions %>% filter(diag_hosp == "R0A07"),
             aes(xintercept = last_icon + 0.5),
             linetype = "dashed", colour = "grey40") +
  geom_line() +
  geom_point(size = 3) +
  scale_colour_manual(values = c(icon = "steelblue", sarah_rapid = "coral")) +
  labs(title = "Colon cancer waiting times", x = NULL, y = "Mean dx_dtt (days)") +
  theme_minimal()



# ── Plot: compare any two years across all hospitals ──────────────────────
year1 <- 2020
year2 <- 2021

compare_data %>%
  filter(ydiag %in% c(year1, year2)) %>%
  select(diag_hosp, ydiag, mean_dx_dtt, n) %>%
  pivot_wider(names_from = ydiag, values_from = c(mean_dx_dtt, n),
              names_sep = "_") %>%
  rename(wait_y1 = paste0("mean_dx_dtt_", year1),
         wait_y2 = paste0("mean_dx_dtt_", year2),
         n_y1    = paste0("n_", year1),
         n_y2    = paste0("n_", year2)) %>%
  filter(!is.na(wait_y1) & !is.na(wait_y2)) %>%
  ggplot(aes(x = wait_y1, y = wait_y2)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(), alpha = 0.75) +
  labs(
    title    = glue::glue("Mean colon cancer wait: {year1} vs {year2}"),
    subtitle = "Above dashed line = longer waits in later year | point size = case volume",
    x        = glue::glue("Mean dx_dtt {year1} (days)"),
    y        = glue::glue("Mean dx_dtt {year2} (days)"),
  ) +
  theme_minimal(base_size = 12)+
  coord_fixed(xlim = c(0, 80), ylim = c(0, 80))



################################################

# Quantify how many cases each filter removes
ncras_hes %>%
  filter(!is.na(tx_date)) %>%
  left_join(cwt_colon_surgery, by = "pseudo_patientid") %>%
  mutate(tx_date_diff = abs(as.numeric(tx_date - cwt_tx_date))) %>%
  summarise(
    n_total          = n(),
    n_no_cwt_match   = sum(is.na(tx_date_diff)),
    n_outside_5days  = sum(!is.na(tx_date_diff) & tx_date_diff > 5),
    n_dtt_zero       = sum(wt_dx_to_dtt == 0, na.rm = TRUE),
    n_dtt_negative   = sum(wt_dx_to_dtt < 0,  na.rm = TRUE),
    n_seq_fail       = sum(!seq_ok, na.rm = TRUE)
  )

test_new <- compare_data %>%
  filter(ydiag %in% c(2018,2019)) %>%
  select(diag_hosp, ydiag, mean_dx_dtt, n) %>%
  filter(!is.na(mean_dx_dtt)) %>%
  group_by(diag_hosp) %>%
  filter(n()==2) %>%
  ungroup() %>%
  pivot_wider(names_from ="ydiag", values_from = c("mean_dx_dtt","n"))

plot(test_new$mean_dx_dtt_2018, test_new$mean_dx_dtt_2019)

test_new %>% View()
cor(test_new$mean_dx_dtt_2018, test_new$mean_dx_dtt_2019)

