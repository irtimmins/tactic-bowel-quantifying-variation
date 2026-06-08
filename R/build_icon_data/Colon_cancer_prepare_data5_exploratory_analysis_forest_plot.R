#############################


# Mean waiting times by covariate group + difference from reference
# Forest plot: 3 parallel panels (wt_dx_to_tx, wt_dx_to_dtt, wt_dtt_to_tx)

library(tidyverse)
library(ggplot2)
library(patchwork)

# ?????? Reference categories ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# Sex:   1 (Male)
# Age:   60-69
# IMD:   Q3
# Route: TWW

# ?????? 1. Build mean + 95% CI by group for each outcome ???????????????????????????????????????????????????????????????????????????

wt_vars   <- c("wt_dx_to_tx", "wt_dx_to_dtt", "wt_dtt_to_tx")
wt_labels <- c("Dx ??? Treatment", "Dx ??? DTT", "DTT ??? Treatment")

df2 <- colon_cohort %>%
  filter(
    !is.na(wt_dx_to_tx), !is.na(wt_dx_to_dtt), !is.na(wt_dtt_to_tx),
    wt_dx_to_tx >= 0, wt_dx_to_dtt >= 0, wt_dtt_to_tx >= 0
  ) %>%
  mutate(
    sex       = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    age_group = cut(agediag, breaks = c(0, 49, 59, 69, 79, Inf),
                    labels = c("<50", "50-59", "60-69", "70-79", "80+"), right = TRUE),
    imd_q     = factor(NHSE_reversed_imd_quintile_lsoas,
                       levels = c("1 - most deprived", "2", "3", "4", "5 - least deprived"),
                       labels = c("Q1 (most deprived)", "Q2", "Q3", "Q4", "Q5 (least deprived)")),
    route     = factor(route_combined)
  )

# Function: mean + 95% CI for one grouping variable, then subtract reference mean
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
      ref_mean = mean[group == ref_level],
      diff     = mean - ref_mean,
      lo       = diff - 1.96 * se,   # approximate; ignores uncertainty in ref
      hi       = diff + 1.96 * se,
      is_ref   = group == ref_level,
      covariate = group_col,
      outcome   = outcome
    )
}

covariates <- list(
  list(col = "sex",       ref = "Male"),
  list(col = "age_group", ref = "60-69"),
  list(col = "imd_q",     ref = "Q3"),
  list(col = "route",     ref = "TWW")
)

results <- map_dfr(wt_vars, function(wt) {
  map_dfr(covariates, function(cv) {
    mean_diff(df2, cv$col, cv$ref, wt)
  })
}) %>%
  mutate(outcome = factor(outcome, levels = wt_vars, labels = wt_labels))

# ?????? 2. Ordering: covariate blocks with a blank separator row ??????????????????????????????????????????????????????

block_order <- c(
  "Male", "Female",
  "<50", "50-59", "60-69", "70-79", "80+",
  "Q1 (most deprived)", "Q2", "Q3", "Q4", "Q5 (least deprived)",
  "GP referral", "Inpatient elective", "Other outpatient", "Screening", "TWW"
)

# Section header labels (inserted as blank rows for spacing in the y-axis)
section_labels <- c(
  "Male"               = "?????? Sex ??????",
  "<50"                = "?????? Age group ??????",
  "Q1 (most deprived)" = "?????? IMD quintile ??????",
  "GP referral"        = "?????? Route to diagnosis ??????"
)

results <- results %>%
  filter(group %in% block_order) %>%
  mutate(group = factor(group, levels = rev(block_order)))

# ?????? 3. Plot ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# After running results, get shared x limits across all three outcomes
x_limits <- c(
  min(results$lo,   na.rm = TRUE),
  max(results$hi,   na.rm = TRUE)
)
# Add a bit of room for the text labels on the right
x_limits[2] <- x_limits[2] * 1.25

make_panel <- function(data, outcome_label, show_y = TRUE) {
  d <- data %>% filter(outcome == outcome_label)
  
  ggplot(d, aes(x = diff, y = group, colour = is_ref)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.6) +
    geom_point(aes(size = is_ref), shape = 18) +
    geom_text(aes(label = ifelse(is_ref, paste0(round(mean), "d (ref)"),
                                 paste0(ifelse(diff >= 0, "+", ""), round(diff), "d"))),
              hjust = -0.15, size = 2.7, colour = "grey20") +
    scale_colour_manual(values = c("TRUE" = "grey60", "FALSE" = "#2166ac")) +
    scale_size_manual(values = c("TRUE" = 2, "FALSE" = 3)) +
    scale_x_continuous(limits = x_limits) +      # <-- shared axis
    facet_grid(covariate ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(covariate = c(
                 sex       = "Sex",
                 age_group = "Age group",
                 imd_q     = "IMD quintile",
                 route     = "Route to diagnosis"
               ))) +
    labs(title = outcome_label, x = "Difference in mean days from reference", y = NULL) +
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

p1 <- make_panel(results, "Dx ??? Treatment",  show_y = TRUE)
p2 <- make_panel(results, "Dx ??? DTT",        show_y = FALSE)
p3 <- make_panel(results, "DTT ??? Treatment", show_y = FALSE)

combined <- p1 + p2 + p3 +
  plot_annotation(
    title    = "Mean waiting times: unadjusted difference from reference category",
    subtitle = "Reference: Male | Age 60-69 | IMD Q3 | TWW route  .  Error bars = ±1.96 SE",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

print(combined)
