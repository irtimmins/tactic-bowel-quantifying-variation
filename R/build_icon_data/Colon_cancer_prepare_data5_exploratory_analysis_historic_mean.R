# =============================================================================
# Colon Cancer Surgery Waiting Times
# Empirical Bayes hospital-level analysis of Dx→DTT
# Period 1: 2015-2019 | Period 2: 2020-2022
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(lme4)
library(lubridate)
library(knitr)

# ── 0. LOAD AND PREPARE ───────────────────────────────────────────────────────

colon_cohort <- readRDS(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_cci_2015_2022.rds"
)

df <- colon_cohort %>%
  filter(
    !is.na(wt_dx_to_tx), !is.na(wt_dx_to_dtt), !is.na(wt_dtt_to_tx),
    wt_dx_to_tx  > 0, wt_dx_to_tx <= 180,
    wt_dx_to_dtt >= 0,
    wt_dtt_to_tx >= 0
  ) %>%
  mutate(
    stage         = factor(stage, levels = 1:3,
                           labels = c("Stage I", "Stage II", "Stage III")),
    route_combined = factor(route_combined),
    sex            = factor(sex, labels = c("Male", "Female")),
    imd_q          = factor(NHSE_reversed_imd_quintile_lsoas,
                            levels = c("1 - most deprived","2","3","4","5 - least deprived"),
                            labels = c("Q1 Most\ndeprived","Q2","Q3","Q4","Q5 Least\ndeprived")),
    ydiag          = as.integer(ydiag),
    age_group      = cut(agediag, breaks = c(0,49,59,69,79,Inf),
                         labels = c("<50","50-59","60-69","70-79","80+"),
                         right  = TRUE)
  )
#

# ── 1. EB DATASET ─────────────────────────────────────────────────────────────

df_eb <- df %>%
  filter(!is.na(wt_dx_to_dtt), !is.na(diag_hosp)) %>%
  mutate(
    season    = case_when(
      month(diagmdy) %in% c(12,1,2)  ~ "Winter",
      month(diagmdy) %in% c(3,4,5)   ~ "Spring",
      month(diagmdy) %in% c(6,7,8)   ~ "Summer",
      month(diagmdy) %in% c(9,10,11) ~ "Autumn"
    ),
    season    = factor(season, levels = c("Spring","Summer","Autumn","Winter")),
    ydiag     = factor(ydiag),
    diag_hosp = factor(diag_hosp)
  )

# Hospital inclusion: present in all 8 years with >=10 cases per year
hosp_min <- df_eb %>%
  count(diag_hosp, ydiag) %>%
  group_by(diag_hosp) %>%
  summarise(
    n_years    = n_distinct(ydiag),
    min_annual = min(n),
    .groups    = "drop"
  ) %>%
  filter(n_years == 8, min_annual >= 10) %>%
  pull(diag_hosp)

cat("Hospitals meeting inclusion criteria:", length(hosp_min), "\n")

df_eb <- df_eb %>%
  filter(diag_hosp %in% hosp_min) %>%
  mutate(diag_hosp = droplevels(diag_hosp))


# ── 2. PERIOD 1: 2015-2019 ────────────────────────────────────────────────────

df_p1 <- df_eb %>%
  filter(ydiag %in% as.character(2015:2020)) %>%
  mutate(ydiag = droplevels(ydiag))

cat("Period 1: n =", nrow(df_p1), "| Hospitals =", n_distinct(df_p1$diag_hosp), "\n")

m1 <- lmer(
  wt_dx_to_dtt ~ age_group + cci_group + ydiag + season + (1 | diag_hosp),
  data    = df_p1,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

# Extract EB estimates
re_obj_p1    <- ranef(m1, condVar = TRUE)$diag_hosp
grand_mean_p1 <- fixef(m1)["(Intercept)"]

eb_p1 <- tibble(
  diag_hosp  = rownames(re_obj_p1),
  re_p1      = re_obj_p1[, 1],
  re_p1_se   = sqrt(as.numeric(attr(re_obj_p1, "postVar")[1, 1, ]))
) %>%
  mutate(
    eb_mean_p1 = grand_mean_p1 + re_p1,
    eb_lo_p1   = grand_mean_p1 + re_p1 - 1.96 * re_p1_se,
    eb_hi_p1   = grand_mean_p1 + re_p1 + 1.96 * re_p1_se,
    n_p1       = as.integer(table(df_p1$diag_hosp)[diag_hosp])
  )

cat("Period 1 EB mean range:", round(range(eb_p1$eb_mean_p1), 1), "\n")


# ── 3. PERIOD 2: 2020-2022 ────────────────────────────────────────────────────

df_p2 <- df_eb %>%
  filter(ydiag %in% as.character(2021:2022)) %>%
  mutate(ydiag = droplevels(factor(ydiag))) %>%
  left_join(eb_p1 %>% select(diag_hosp, eb_mean_p1), by = "diag_hosp")

cat("Period 2: n =", nrow(df_p2), "| Hospitals =", n_distinct(df_p2$diag_hosp), "\n")


# ── 4. RAW MEANS CATERPILLAR (Period 2) ───────────────────────────────────────

raw_means <- df_p2 %>%
  group_by(diag_hosp) %>%
  summarise(
    raw_mean = mean(wt_dx_to_dtt, na.rm = TRUE),
    raw_se   = sd(wt_dx_to_dtt,   na.rm = TRUE) / sqrt(n()),
    n        = n(),
    .groups  = "drop"
  ) %>%
  mutate(
    raw_lo    = raw_mean - 1.96 * raw_se,
    raw_hi    = raw_mean + 1.96 * raw_se,
    grand_raw = mean(raw_mean),
    diag_hosp = fct_reorder(diag_hosp, raw_mean),
    sig       = case_when(
      raw_lo > grand_raw ~ "Above average",
      raw_hi < grand_raw ~ "Below average",
      TRUE               ~ "Not significant"
    ),
    sig = factor(sig, levels = c("Above average","Not significant","Below average"))
  )

ggplot(raw_means, aes(x = raw_mean, y = diag_hosp, colour = sig)) +
  geom_vline(xintercept = mean(raw_means$raw_mean), linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = raw_lo, xmax = raw_hi),
                 height = 0, linewidth = 0.4, alpha = 0.6) +
  geom_point(aes(size = n), alpha = 0.8) +
  scale_colour_manual(values = c("Above average" = "#E15759",
                                 "Not significant" = "grey60",
                                 "Below average"   = "#4E79A7")) +
  scale_size_continuous(range = c(1,5), name = "N (2020-22)") +
  labs(title    = "Raw mean Dx\u2192DTT per hospital (2020-22)",
       subtitle = "Unadjusted; ranked by mean; error bars = 95% CI; dashed = grand mean",
       x = "Mean days (Dx \u2192 DTT)", y = NULL, colour = NULL) +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        axis.text.y = element_text(size = 6),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")


# ── 5. MODEL A: Period 2 without historic adjustment ─────────────────────────

m2a <- lmer(
  wt_dx_to_dtt ~ age_group + cci_group + ydiag + season + (1 | diag_hosp),
  data    = df_p2,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

re_obj_a     <- ranef(m2a, condVar = TRUE)$diag_hosp
grand_mean_a <- fixef(m2a)["(Intercept)"]

eb_a <- tibble(
  diag_hosp = rownames(re_obj_a),
  re_a      = re_obj_a[, 1],
  re_a_se   = sqrt(as.numeric(attr(re_obj_a, "postVar")[1, 1, ]))
) %>%
  mutate(
    eb_mean_a = re_a,
    eb_lo_a   = re_a - 1.96 * re_a_se,
    eb_hi_a   = re_a + 1.96 * re_a_se
  )

cat("Model A — SD of EB means:", round(sd(eb_a$eb_mean_a), 2), "\n")
cat("Model A — ICC:", round(
  as.data.frame(VarCorr(m2a))$vcov[1] / sum(as.data.frame(VarCorr(m2a))$vcov), 3), "\n")


# ── 6. MODEL B: Period 2 with historic adjustment ─────────────────────────────

m2b <- lmer(
  wt_dx_to_dtt ~ age_group + cci_group + ydiag + season + eb_mean_p1 +
    (1 | diag_hosp),
  data    = df_p2,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

re_obj_b     <- ranef(m2b, condVar = TRUE)$diag_hosp
grand_mean_b <- fixef(m2b)["(Intercept)"]

eb_b <- tibble(
  diag_hosp = rownames(re_obj_b),
  re_b      = re_obj_b[, 1],
  re_b_se   = sqrt(as.numeric(attr(re_obj_b, "postVar")[1, 1, ]))
) %>%
  mutate(
    eb_mean_b = re_b,
    eb_lo_b   = re_b - 1.96 * re_b_se,
    eb_hi_b   = re_b + 1.96 * re_b_se
  )

cat("Model B — SD of EB means:", round(sd(eb_b$eb_mean_b), 2), "\n")


# ── 7. COMPARE MODEL A vs MODEL B ────────────────────────────────────────────

eb_compare <- eb_a %>%
  select(diag_hosp, eb_mean_a, eb_lo_a, eb_hi_a) %>%
  inner_join(eb_b %>% select(diag_hosp, eb_mean_b, eb_lo_b, eb_hi_b),
             by = "diag_hosp") %>%
  inner_join(eb_p1 %>% select(diag_hosp, eb_mean_p1, n_p1),
             by = "diag_hosp") %>%
  left_join(df_p2 %>% count(diag_hosp) %>% rename(n_p2 = n),
            by = "diag_hosp") %>%
  mutate(
    shift     = eb_mean_b - eb_mean_a,
    shift_dir = case_when(
      shift >  2 ~ "Higher after baseline adjustment",
      shift < -2 ~ "Lower after baseline adjustment",
      TRUE       ~ "Minimal change"
    ),
    shift_dir = factor(shift_dir, levels = c(
      "Higher after baseline adjustment",
      "Lower after baseline adjustment",
      "Minimal change"
    ))
  )

cat("Hospitals in comparison:", nrow(eb_compare), "\n")
cat("Pearson r (A vs B):",
    round(cor(eb_compare$eb_mean_a, eb_compare$eb_mean_b), 3), "\n")
cat("Spearman r (A vs B):",
    round(cor(eb_compare$eb_mean_a, eb_compare$eb_mean_b, method = "spearman"), 3), "\n")



# =============================================================================
# Model A and B variance decomposition
# =============================================================================

# Variance components
vc_a <- as.data.frame(VarCorr(m2a))
vc_b <- as.data.frame(VarCorr(m2b))

# Hospital and residual variance
var_hosp_a  <- vc_a$vcov[1]
var_resid_a <- vc_a$vcov[2]
var_total_a <- var_hosp_a + var_resid_a

var_hosp_b  <- vc_b$vcov[1]
var_resid_b <- vc_b$vcov[2]
var_total_b <- var_hosp_b + var_resid_b

cat("=== MODEL A: Period 2, no historic adjustment ===\n")
cat("Hospital-level variance:      ", round(var_hosp_a,  2), "\n")
cat("Hospital-level SD:            ", round(sqrt(var_hosp_a), 2), "\n")
cat("Residual variance:            ", round(var_resid_a, 2), "\n")
cat("Total variance:               ", round(var_total_a, 2), "\n")
cat("ICC (% variance at hospital): ", round(100 * var_hosp_a / var_total_a, 1), "%\n")

cat("\n=== MODEL B: Period 2, EB historic adjustment ===\n")
cat("Hospital-level variance:      ", round(var_hosp_b,  2), "\n")
cat("Hospital-level SD:            ", round(sqrt(var_hosp_b), 2), "\n")
cat("Residual variance:            ", round(var_resid_b, 2), "\n")
cat("Total variance:               ", round(var_total_b, 2), "\n")
cat("ICC (% variance at hospital): ", round(100 * var_hosp_b / var_total_b, 1), "%\n")

cat("\n=== VARIANCE EXPLAINED BY HISTORIC EB MEAN ===\n")
cat("Reduction in hospital variance (B vs A): ",
    round(100 * (var_hosp_a - var_hosp_b) / var_hosp_a, 1), "%\n")
cat("Reduction in total variance (B vs A):    ",
    round(100 * (var_total_a - var_total_b) / var_total_a, 1), "%\n")

cat("\n=== FIXED EFFECTS R² (Nakagawa & Schielzeth) ===\n")
# Marginal R2 = variance explained by fixed effects only
# Conditional R2 = variance explained by fixed + random effects
# Using performance package if available, otherwise manual

if (requireNamespace("performance", quietly = TRUE)) {
  library(performance)
  cat("\nModel A:\n"); print(r2(m2a))
  cat("\nModel B:\n"); print(r2(m2b))
} else {
  # Manual Nakagawa R2
  # Marginal: var(fixed) / (var(fixed) + var(random) + var(residual))
  # Conditional: (var(fixed) + var(random)) / (var(fixed) + var(random) + var(residual))
  
  var_fixed_a <- var(predict(m2a, re.form = NA))
  var_fixed_b <- var(predict(m2b, re.form = NA))
  
  r2_marginal_a    <- var_fixed_a / (var_fixed_a + var_hosp_a + var_resid_a)
  r2_conditional_a <- (var_fixed_a + var_hosp_a) / (var_fixed_a + var_hosp_a + var_resid_a)
  
  r2_marginal_b    <- var_fixed_b / (var_fixed_b + var_hosp_b + var_resid_b)
  r2_conditional_b <- (var_fixed_b + var_hosp_b) / (var_fixed_b + var_hosp_b + var_resid_b)
  
  cat("\nModel A:\n")
  cat("  Marginal R2 (fixed effects only):        ", round(r2_marginal_a,    3), "\n")
  cat("  Conditional R2 (fixed + random effects): ", round(r2_conditional_a, 3), "\n")
  
  cat("\nModel B:\n")
  cat("  Marginal R2 (fixed effects only):        ", round(r2_marginal_b,    3), "\n")
  cat("  Conditional R2 (fixed + random effects): ", round(r2_conditional_b, 3), "\n")
  
  cat("\nIncrease in marginal R2 due to historic EB mean (B - A):",
      round(r2_marginal_b - r2_marginal_a, 3), "\n")
}


# ── 8. SCATTER: Model A vs Model B ───────────────────────────────────────────

ax_min <- min(c(eb_compare$eb_mean_a, eb_compare$eb_mean_b)) * 0.90
ax_max <- max(c(eb_compare$eb_mean_a, eb_compare$eb_mean_b)) * 1.11

ggplot(eb_compare, aes(x = eb_mean_a, y = eb_mean_b)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              colour = "grey50", linewidth = 0.5) +
  geom_errorbar(aes(ymin = eb_lo_b, ymax = eb_hi_b),
                width = 0, linewidth = 0.3, alpha = 0.2) +
  geom_errorbarh(aes(xmin = eb_lo_a, xmax = eb_hi_a),
                 height = 0, linewidth = 0.3, alpha = 0.2) +
  geom_point(alpha = 0.8, colour = "grey30", size = 3) +
  coord_equal(xlim = c(ax_min, ax_max), ylim = c(ax_min, ax_max)) +
  labs(
    title    = "Hospital EB mean Dx\u2192DTT (2020-22): Model A vs Model B",
    subtitle = "Model A: no historic adjustment | Model B: adjusted for 2015-19 EB mean\nAbove diagonal = longer wait after adjustment",
    x        = "Model A: EB mean days (no historic adjustment)",
    y        = "Model B: EB mean days (adjusted for historic performance)",
    colour   = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        legend.position = "bottom")


# ── 9. RANK CHANGE TABLE ──────────────────────────────────────────────────────

eb_compare %>%
  mutate(
    rank_a      = rank(eb_mean_a),
    rank_b      = rank(eb_mean_b),
    rank_change = rank_b - rank_a
  ) %>%
  arrange(desc(abs(rank_change))) %>%
  select(diag_hosp, eb_mean_p1, eb_mean_a, eb_mean_b,
         shift, rank_a, rank_b, rank_change, n_p1, n_p2) %>%
  print(n = 20)


# Is the historic mean coefficient in Model C significant and in expected direction?
summary(m2c)

# Specifically the raw_mean_p1 fixed effect
fixef(m2c)["raw_mean_p1"]
confint(m2c, parm = "raw_mean_p1", method = "Wald")




# ── 10. raw meam

# Add raw historic mean per hospital from Period 1
raw_mean_p1 <- df_p1 %>%
  group_by(diag_hosp) %>%
  summarise(raw_mean_p1 = mean(wt_dx_to_dtt, na.rm = TRUE),
            .groups = "drop")

# Join to df_p2
df_p2_raw <- df_p2 %>%
  left_join(raw_mean_p1, by = "diag_hosp")

# Model C: Period 2 adjusted for raw historic mean instead of EB mean
m2c <- lmer(
  wt_dx_to_dtt ~ age_group + cci_group + ydiag + season + raw_mean_p1 +
    (1 | diag_hosp),
  data    = df_p2_raw,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

re_obj_c     <- ranef(m2c, condVar = TRUE)$diag_hosp
grand_mean_c <- fixef(m2c)["(Intercept)"]

eb_c <- tibble(
  diag_hosp = rownames(re_obj_c),
  re_c      = re_obj_c[, 1],
  re_c_se   = sqrt(as.numeric(attr(re_obj_c, "postVar")[1, 1, ]))
) %>%
  mutate(
    eb_mean_c = grand_mean_c + re_c,
    eb_lo_c   = grand_mean_c + re_c - 1.96 * re_c_se,
    eb_hi_c   = grand_mean_c + re_c + 1.96 * re_c_se
  )

cat("Model C (raw historic) — SD of EB means:", round(sd(eb_c$eb_mean_c), 2), "\n")

# Compare all three models
eb_all <- eb_compare %>%
  inner_join(eb_c %>% select(diag_hosp, eb_mean_c, eb_lo_c, eb_hi_c),
             by = "diag_hosp") %>%
  inner_join(raw_mean_p1, by = "diag_hosp")

cat("\nCorrelations:\n")
cat("A vs B (EB historic):",
    round(cor(eb_all$eb_mean_a, eb_all$eb_mean_b), 3), "\n")
cat("A vs C (raw historic):",
    round(cor(eb_all$eb_mean_a, eb_all$eb_mean_c), 3), "\n")
cat("B vs C (EB vs raw):",
    round(cor(eb_all$eb_mean_b, eb_all$eb_mean_c), 3), "\n")

# Scatter: B vs C — do EB and raw historic adjustments give similar results?
ax_min <- min(c(eb_all$eb_mean_b, eb_all$eb_mean_c)) * 0.95
ax_max <- max(c(eb_all$eb_mean_b, eb_all$eb_mean_c)) * 1.05

ggplot(eb_all, aes(x = eb_mean_b, y = eb_mean_c)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              colour = "grey50", linewidth = 0.5) +
  geom_point(aes(size = n_p2), alpha = 0.7, colour = "#4E79A7") +
  coord_equal(xlim = c(ax_min, ax_max), ylim = c(ax_min, ax_max)) +
  labs(
    title    = "Model B (EB historic) vs Model C (raw historic mean)",
    subtitle = "Comparing two approaches to adjusting for Period 1 baseline performance",
    x        = "Model B: EB mean days (EB historic adjustment)",
    y        = "Model C: EB mean days (raw historic adjustment)",
    size     = "N (2020-22)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        legend.position = "bottom")






######################################################


# ── Raw annual means per hospital 2015-2019 ───────────────────────────────────

raw_annual_p1 <- df_p1 %>%
  group_by(diag_hosp, ydiag) %>%
  summarise(mean_wt = mean(wt_dx_to_dtt, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(ydiag = paste0("mean_", ydiag)) %>%
  pivot_wider(names_from = ydiag, values_from = mean_wt)

# Check
print(raw_annual_p1, n = 5)

# Join to df_p2
df_p2_annual <- df_p2 %>%
  left_join(raw_annual_p1, by = "diag_hosp")

# Check for NAs — should be none given hosp_min criterion
df_p2_annual %>%
  summarise(across(starts_with("mean_"), ~sum(is.na(.))))

# ── Model D: Period 2 with all 5 annual historic means ───────────────────────

m2d <- lmer(
  wt_dx_to_dtt ~ age_group + cci_group + ydiag + season +
    mean_2015 + mean_2016 + mean_2017 + mean_2018 + mean_2019 +
    (1 | diag_hosp),
  data    = df_p2_annual,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(m2d)

# Variance decomposition
vc_d        <- as.data.frame(VarCorr(m2d))
var_hosp_d  <- vc_d$vcov[1]
var_resid_d <- vc_d$vcov[2]
var_total_d <- var_hosp_d + var_resid_d

cat("=== MODEL D: Period 2, all 5 annual historic means ===\n")
cat("Hospital-level variance:      ", round(var_hosp_d,  2), "\n")
cat("Hospital-level SD:            ", round(sqrt(var_hosp_d), 2), "\n")
cat("Residual variance:            ", round(var_resid_d, 2), "\n")
cat("ICC (% variance at hospital): ", round(100 * var_hosp_d / var_total_d, 1), "%\n")
cat("Reduction in hospital variance vs Model A:",
    round(100 * (var_hosp_a - var_hosp_d) / var_hosp_a, 1), "%\n")
cat("Reduction in hospital variance vs Model B:",
    round(100 * (var_hosp_b - var_hosp_d) / var_hosp_b, 1), "%\n")

# Fixed effect coefficients on historic years
cat("\nHistoric year coefficients:\n")
fixef(m2d)[grep("mean_", names(fixef(m2d)))]

# Extract EB estimates
re_obj_d     <- ranef(m2d, condVar = TRUE)$diag_hosp
grand_mean_d <- fixef(m2d)["(Intercept)"]

eb_d <- tibble(
  diag_hosp = rownames(re_obj_d),
  re_d      = re_obj_d[, 1],
  re_d_se   = sqrt(as.numeric(attr(re_obj_d, "postVar")[1, 1, ]))
) %>%
  mutate(
    eb_mean_d = grand_mean_d + re_d,
    eb_lo_d   = grand_mean_d + re_d - 1.96 * re_d_se,
    eb_hi_d   = grand_mean_d + re_d + 1.96 * re_d_se
  )

cat("\nModel D — SD of EB means:", round(sd(eb_d$eb_mean_d), 2), "\n")

# Compare A, B, D
eb_all3 <- eb_compare %>%
  inner_join(eb_d %>% select(diag_hosp, eb_mean_d), by = "diag_hosp")

cat("\nCorrelations:\n")
cat("A vs B:", round(cor(eb_all3$eb_mean_a, eb_all3$eb_mean_b), 3), "\n")
cat("A vs D:", round(cor(eb_all3$eb_mean_a, eb_all3$eb_mean_d), 3), "\n")
cat("B vs D:", round(cor(eb_all3$eb_mean_b, eb_all3$eb_mean_d), 3), "\n")

##################################################################



# =============================================================================
# Posterior probability of being in top 5%, 10%, 20%, 50% of providers
# Based on Model A and Model B EB posterior distributions
# Lower EB mean = shorter wait = better performance
# =============================================================================

library(tidyverse)

# Number of hospitals
n_hosp <- nrow(eb_a)

# Thresholds (top X% = bottom X% of waiting time distribution)
top_pct <- c(0.05, 0.10, 0.20, 0.50)

# Function: probability hospital i is in top X% 
# Uses Monte Carlo simulation from posterior distributions
set.seed(42)
n_sim <- 10000

simulate_ranks <- function(eb_means, eb_ses, n_sim = 10000) {
  n_hosp <- length(eb_means)
  # Draw n_sim samples from each hospital's posterior
  draws <- mapply(function(mu, se) rnorm(n_sim, mu, se),
                  eb_means, eb_ses,
                  SIMPLIFY = TRUE)  # n_sim x n_hosp matrix
  # For each simulation, rank hospitals (rank 1 = shortest wait = best)
  ranks <- apply(draws, 1, function(x) rank(x))  # n_hosp x n_sim
  return(ranks)
}

# ── Model A ───────────────────────────────────────────────────────────────────

ranks_a <- simulate_ranks(eb_a$eb_mean_a, eb_a$re_a_se)

prob_top_a <- map_dfr(top_pct, function(p) {
  threshold <- ceiling(n_hosp * p)
  probs <- rowMeans(ranks_a <= threshold)
  tibble(
    diag_hosp = eb_a$diag_hosp,
    threshold = paste0("Top ", p * 100, "%"),
    prob      = probs
  )
}) %>%
  pivot_wider(names_from = threshold, values_from = prob) %>%
  left_join(eb_a %>% select(diag_hosp, eb_mean_a), by = "diag_hosp") %>%
  arrange(eb_mean_a)

cat("=== Model A: Posterior probabilities of top performance ===\n")
print(prob_top_a, n = Inf)


# ── Model B ───────────────────────────────────────────────────────────────────

ranks_b <- simulate_ranks(eb_b$eb_mean_b, eb_b$re_b_se)

prob_top_b <- map_dfr(top_pct, function(p) {
  threshold <- ceiling(n_hosp * p)
  probs <- rowMeans(ranks_b <= threshold)
  tibble(
    diag_hosp = eb_b$diag_hosp,
    threshold = paste0("Top ", p * 100, "%"),
    prob      = probs
  )
}) %>%
  pivot_wider(names_from = threshold, values_from = prob) %>%
  left_join(eb_b %>% select(diag_hosp, eb_mean_b), by = "diag_hosp") %>%
  arrange(eb_mean_b)

cat("\n=== Model B: Posterior probabilities of top performance ===\n")
print(prob_top_b, n = 20)


# ── COMBINED TABLE ────────────────────────────────────────────────────────────

prob_combined <- prob_top_a %>%
  select(diag_hosp, eb_mean_a,
         `A: Top 5%`  = `Top 5%`,
         `A: Top 10%` = `Top 10%`,
         `A: Top 20%` = `Top 20%`,
         `A: Top 50%` = `Top 50%`) %>%
  inner_join(
    prob_top_b %>%
      select(diag_hosp, eb_mean_b,
             `B: Top 5%`  = `Top 5%`,
             `B: Top 10%` = `Top 10%`,
             `B: Top 20%` = `Top 20%`,
             `B: Top 50%` = `Top 50%`),
    by = "diag_hosp"
  ) %>%
  arrange(eb_mean_a)

print(prob_combined, n = Inf)



compare_rankings <- prob_combined %>%
  mutate(
    rank_a      = rank(-`A: Top 20%`, ties.method = "min"),
    rank_b      = rank(-`B: Top 20%`, ties.method = "min"),
    rank_change = rank_a - rank_b
  ) %>%
  arrange(rank_b) %>%
  select(
    Hospital        = diag_hosp,
    `EB mean (B)`   = eb_mean_b,
    `EB mean (A)`   = eb_mean_a,
    `P(Top 10%) A`  = `A: Top 10%`,
    `P(Top 10%) B`  = `B: Top 10%`,
    `P(Top 20%) A`  = `A: Top 20%`,
    `P(Top 20%) B`  = `B: Top 20%`,
    `Rank (B)`      = rank_b,
    `Rank (A)`      = rank_a,
    `Rank change`   = rank_change
  ) %>%
  mutate(
    `EB mean (B)`  = formatC(`EB mean (B)`,  format = "f", digits = 1),
    `EB mean (A)`  = formatC(`EB mean (A)`,  format = "f", digits = 1),
    `P(Top 10%) A` = formatC(`P(Top 10%) A`, format = "f", digits = 2),
    `P(Top 10%) B` = formatC(`P(Top 10%) B`, format = "f", digits = 2),
    `P(Top 20%) A` = formatC(`P(Top 20%) A`, format = "f", digits = 2),
    `P(Top 20%) B` = formatC(`P(Top 20%) B`, format = "f", digits = 2)
  )

print(compare_rankings, n = Inf)

print(compare_rankings, n = Inf)
print(compare_rankings, n = Inf)
  
#View()

write_csv(compare_rankings, "D:/Projects/#2045_ICON_TACTIC/Project3_MSc_ses/Prepare_data/compare_rankings_A_B.csv")


# ── PLOT: probability of being in top 20% — Model A vs Model B ───────────────

prob_combined %>%
  mutate(diag_hosp = fct_reorder(diag_hosp, `A: Top 20%`)) %>%
  ggplot(aes(x = `A: Top 20%`, y = `B: Top 20%`,
             colour = `A: Top 20%` - `B: Top 20%`)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              colour = "grey50", linewidth = 0.5) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_gradient2(low = "#E15759", mid = "grey80", high = "#4E79A7",
                         midpoint = 0,
                         name = "Change in prob\n(A - B)") +
  labs(
    title    = "Probability of being in top 20% of providers",
    subtitle = "Model A (no historic adjustment) vs Model B (EB historic adjustment)\nAbove diagonal = higher probability after adjustment",
    x        = "Model A: P(top 20%)",
    y        = "Model B: P(top 20%)"
  ) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        legend.position = "right")


# ── CATERPILLAR: P(top 20%) ranked by Model A ────────────────────────────────

prob_combined %>%
  pivot_longer(cols = c(`A: Top 20%`, `B: Top 20%`),
               names_to = "model", values_to = "prob") %>%
  mutate(diag_hosp = fct_reorder(diag_hosp, 
                                 if_else(model == "A: Top 20%", prob, NA_real_),
                                 .fun = function(x) mean(x, na.rm = TRUE))) %>%
  ggplot(aes(x = prob, y = diag_hosp, colour = model)) +
  geom_vline(xintercept = 0.20, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_line(aes(group = diag_hosp), colour = "grey80", linewidth = 0.3) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_colour_manual(values = c("A: Top 20%" = "#4E79A7",
                                 "B: Top 20%" = "#E15759")) +
  labs(
    title    = "P(top 20%) per hospital: Model A vs Model B",
    subtitle = "Ranked by Model A probability; dashed = 20% reference",
    x        = "Posterior probability of being in top 20%",
    y        = NULL,
    colour   = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.title         = element_text(face = "bold", size = 12),
        plot.subtitle      = element_text(size = 9, colour = "grey40"),
        axis.text.y        = element_text(size = 6),
        panel.grid.major.y = element_blank(),
        legend.position    = "bottom")

#install.packages("ggvenn")
library(ggvenn)

# Hospitals clearing 80% probability threshold at each band — Model A vs B

venn_data_5 <- list(
  "Model A\n(Top 5%)"  = prob_combined %>% filter(`A: Top 5%`  > 0.8) %>% pull(diag_hosp) %>% as.character(),
  "Model B\n(Top 5%)"  = prob_combined %>% filter(`B: Top 5%`  > 0.8) %>% pull(diag_hosp) %>% as.character()
)

venn_data_10 <- list(
  "Model A\n(Top 10%)" = prob_combined %>% filter(`A: Top 10%` > 0.8) %>% pull(diag_hosp) %>% as.character(),
  "Model B\n(Top 10%)" = prob_combined %>% filter(`B: Top 10%` > 0.8) %>% pull(diag_hosp) %>% as.character()
)

venn_data_20 <- list(
  "Model A\n(Top 20%)" = prob_combined %>% filter(`A: Top 20%` > 0.8) %>% pull(diag_hosp) %>% as.character(),
  "Model B\n(Top 20%)" = prob_combined %>% filter(`B: Top 20%` > 0.8) %>% pull(diag_hosp) %>% as.character()
)

venn_data_50 <- list(
  "Model A\n(Top 50%)" = prob_combined %>% filter(`A: Top 50%` > 0.8) %>% pull(diag_hosp) %>% as.character(),
  "Model B\n(Top 50%)" = prob_combined %>% filter(`B: Top 50%` > 0.8) %>% pull(diag_hosp) %>% as.character()
)

p_v5 <- ggvenn(venn_data_5,
               fill_color  = c("#4E79A7","#E15759"),
               stroke_size = 0.5, set_name_size = 3.5, text_size = 3.5) +
  labs(title = "P(top 5%) > 80%") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

p_v10 <- ggvenn(venn_data_10,
                fill_color  = c("#4E79A7","#E15759"),
                stroke_size = 0.5, set_name_size = 3.5, text_size = 3.5) +
  labs(title = "P(top 10%) > 80%") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

p_v20 <- ggvenn(venn_data_20,
                fill_color  = c("#4E79A7","#E15759"),
                stroke_size = 0.5, set_name_size = 3.5, text_size = 3.5) +
  labs(title = "P(top 20%) > 80%") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

p_v50 <- ggvenn(venn_data_50,
                fill_color  = c("#4E79A7","#E15759"),
                stroke_size = 0.5, set_name_size = 3.5, text_size = 3.5) +
  labs(title = "P(top 50%) > 80%") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

(p_v5 | p_v10) / (p_v20 | p_v50) +
  plot_annotation(
    title    = "Hospitals with >80% posterior probability of top performance",
    subtitle = "Model A: no historic adjustment | Model B: adjusted for 2015-19 EB mean",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

