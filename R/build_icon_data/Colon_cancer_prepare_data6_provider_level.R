
library(readxl)


colon_cohort <- readRDS(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_2015_2022.rds"
)

provider_df <- read_excel("D:/Projects/#2045_ICON_TACTIC/Project3_MSc_ses/Provider_level/NHSHospitals_services_5.3.26.xlsx", 
                                           skip = 1)

length((unique(provider_df$`Hospital Name`)))
length((unique(provider_df$`Hospital site code`)))
#length(unique(colon_cohort))
#View(NHSHospitals_services_5_3_26)
#provider_df

summary(colon_cohort$wt_dx_to_tx)
summary(colon_cohort$wt_dx_to_dtt)

print(provider_df, n = 10, width = Inf)

#View(provider_df)
provider_df %>%
  filter(!is.na(`Hospital site code`)) %>%
  group_by(`Hospital site code`) %>%
  mutate(count = n()) %>%
  arrange(-count, `Hospital site code`) %>%
  print(n = 20)

#colon_cohort$diag_trust

##################################################



summary(provider_df$mean)
hist(provider_df$mean, breaks = 20)

hist(provider_df$`Staff engagement`, breaks = 30)
hist(provider_df$Moral, breaks = 30)
hist(provider_df$`We are a team`, breaks = 30)

###################################################



# =============================================================================
# Provider-level analysis: merge & stratify waiting times
# =============================================================================

library(tidyverse)
library(readxl)

# ?????? 1. COLLAPSE PROVIDER DATA TO TRUST LEVEL ??????????????????????????????????????????????????????????????????????????????????????????????????????
View(provider_df)
provider_trust <- provider_df %>%
  group_by(trust_nacs) %>%
  summarise(
    comprehensive_centre = max(`Comprehensive centre`, na.rm = TRUE),
    teaching_hospital    = max(`Teaching hospitals`,  na.rm = TRUE),
    cqc_rating           = first(na.omit(`Latest Rating`)),
    staff_engagement     = first(na.omit(`Staff engagement`)),
    moral                = first(na.omit(Moral)),
    bed_occ_mean         = first(na.omit(mean)),   # NHS performance = bed occupancy
    .groups = "drop"
  ) %>%
  mutate(
    # Replace -Inf from max() on all-NA columns back to NA
    across(c(comprehensive_centre, teaching_hospital),
           ~ifelse(is.infinite(.), NA, .)),
    # Binary factors
    comprehensive_centre = factor(comprehensive_centre, levels = c(0,1),
                                  labels = c("Non-comprehensive", "Comprehensive")),
    teaching_hospital    = factor(teaching_hospital,    levels = c(0,1),
                                  labels = c("Non-teaching", "Teaching")),
    # CQC ordered factor
    cqc_rating = factor(cqc_rating,
                        levels = c("Inadequate", "Requires Improvement", "Good", "Outstanding"),
                        ordered = TRUE),
    # Staff engagement and morale: trust-level quintiles, Q3 as reference
    staff_eng_cat = cut(staff_engagement,
                        breaks = quantile(staff_engagement, probs = seq(0, 1, 0.2), na.rm = TRUE),
                        labels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"),
                        include.lowest = TRUE),
    moral_cat     = cut(moral,
                        breaks = quantile(moral, probs = seq(0, 1, 0.2), na.rm = TRUE),
                        labels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"),
                        include.lowest = TRUE), # Bed occupancy: cut at 92% and 95%
    bed_occ_cat = case_when(
      is.na(bed_occ_mean)      ~ NA_character_,
      bed_occ_mean >= 0.95     ~ "High (>=95%)",
      TRUE                     ~ "Normal (<95%)"
    ),
    bed_occ_cat = factor(bed_occ_cat,
                         levels = c("Normal (<95%)","High (>=95%)"))
  ) 


# ?????? 2. JOIN TO COHORT ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????

df_prov <- colon_cohort %>%
  filter(
    !is.na(wt_dx_to_tx), !is.na(wt_dx_to_dtt), !is.na(wt_dtt_to_tx),
    wt_dx_to_tx >= 0, wt_dx_to_dtt >= 0, wt_dtt_to_tx >= 0
  ) %>%
  left_join(provider_trust, by = c("diag_trust" = "trust_nacs"))

# Join diagnostics
cat("Cohort n =", nrow(df_prov), "\n")
cat("Matched to provider data:", sum(!is.na(df_prov$cqc_rating)), "\n")
cat("Unmatched trusts:", paste(setdiff(unique(df_prov$diag_trust),
                                       provider_trust$trust_nacs), collapse = ", "), "\n")


# ?????? 3. MEAN DIFFERENCE FUNCTION (reused from earlier) ????????????????????????????????????????????????????????????????????????

wt_vars   <- c("wt_dx_to_tx", "wt_dx_to_dtt", "wt_dtt_to_tx")
wt_labels <- c("Dx ??? Treatment", "Dx ??? DTT", "DTT ??? Treatment")

# Winsorise at 99th percentile
#cap99 <- function(x) pmin(x, quantile(x, 0.99, na.rm = TRUE))
#df_prov <- df_prov %>%
 # mutate(across(all_of(wt_vars), cap99))

mean_diff <- function(data, group_col, ref_level, outcome) {
  data %>%
    filter(!is.na(.data[[group_col]]), !is.na(.data[[outcome]])) %>%
    group_by(group = .data[[group_col]]) %>%
    summarise(
      n    = n(),
      mean = mean(.data[[outcome]], na.rm = TRUE),
      se   = sd(.data[[outcome]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      group    = as.character(group),   # <-- fix: drop factor/ordered distinction
      ref_mean = mean[group == ref_level],
      diff     = mean - ref_mean,
      lo       = diff - 1.96 * se,
      hi       = diff + 1.96 * se,
      is_ref   = group == ref_level,
      covariate = group_col,
      outcome   = outcome
    )
}

# ?????? 4. BUILD RESULTS ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????

covariates <- list(
  list(col = "comprehensive_centre", ref = "Non-comprehensive"),
  list(col = "teaching_hospital",    ref = "Non-teaching"),
  list(col = "cqc_rating",           ref = "Good"),
  list(col = "staff_eng_cat",        ref = "Q3"),
  list(col = "moral_cat",            ref = "Q3"),
  list(col = "bed_occ_cat",          ref = "Normal (<95%)")
)

results_prov <- map_dfr(wt_vars, function(wt) {
  map_dfr(covariates, function(cv) {
    mean_diff(df_prov, cv$col, cv$ref, wt)
  })
}) %>%
  mutate(outcome = factor(outcome, levels = wt_vars, labels = wt_labels))

# Group ordering for y-axis
level_order <- c(
  "Non-comprehensive", "Comprehensive",
  "Non-teaching", "Teaching",
  "Inadequate", "Requires Improvement", "Good", "Outstanding",
  "Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)",  # used by both staff_eng_cat and moral_cat
  "Normal (<95%)", "High (>=95%)"
)

results_prov <- results_prov %>%
  filter(group %in% level_order) %>%
  mutate(group = factor(group, levels = rev(level_order)))


# ?????? 5. FOREST PLOT ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????

x_limits <- c(
  min(results_prov$lo, na.rm = TRUE),
  max(results_prov$hi, na.rm = TRUE) * 1.25
)

make_panel_prov <- function(data, outcome_label, show_y = TRUE) {
  d <- data %>% filter(outcome == outcome_label)
  
  ggplot(d, aes(x = diff, y = group, colour = is_ref)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.6) +
    geom_point(aes(size = is_ref), shape = 18) +
    geom_text(aes(label = ifelse(is_ref,
                                 paste0(round(mean), "d (ref)"),
                                 paste0(ifelse(diff >= 0, "+", ""), round(diff), "d"))),
              hjust = -0.15, size = 2.7, colour = "grey20") +
    scale_colour_manual(values = c("TRUE" = "grey60", "FALSE" = "#2166ac")) +
    scale_size_manual(values   = c("TRUE" = 2,        "FALSE" = 3)) +
    scale_x_continuous(limits  = x_limits) +
    facet_grid(covariate ~ ., scales = "free_y", space = "free_y",
               # in facet labeller:
               labeller = labeller(covariate = c(
                 comprehensive_centre = "Centre type",
                 teaching_hospital    = "Teaching status",
                 cqc_rating           = "CQC rating",
                 staff_eng_cat        = "Staff engagement",
                 moral_cat            = "Morale",
                 bed_occ_cat          = "Bed occupancy"
               ))) +
    labs(title = outcome_label,
         x = "Difference in mean days from reference", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position    = "none",
      strip.text.y       = element_text(angle = 0, hjust = 0, face = "bold", size = 9),
      strip.background   = element_rect(fill = "grey93", colour = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text.y        = if (show_y) element_text(size = 9) else element_blank(),
      axis.ticks.y       = element_blank()
    )
}

p1 <- make_panel_prov(results_prov, "Dx ??? Treatment",  show_y = TRUE)
p2 <- make_panel_prov(results_prov, "Dx ??? DTT",        show_y = FALSE)
p3 <- make_panel_prov(results_prov, "DTT ??? Treatment", show_y = FALSE)

combined_prov <- p1 + p2 + p3 +
  plot_annotation(
    title    = "Waiting times by provider characteristics (unadjusted)",
    subtitle = "Reference: Non-comprehensive | Non-teaching | CQC Good | Staff engagement <90th | Bed occupancy <92%\nError bars = ±1.96 SE",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

print(combined_prov)





























