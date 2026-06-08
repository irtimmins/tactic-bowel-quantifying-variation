# =============================================================================
# Synthetic provider-characteristics workbook for the colon site-level merge
# -----------------------------------------------------------------------------
# Builds a small Excel that mimics NHSHospitals_services_*.xlsx but contains
# only the columns the Stata provider step (06) reads, one row per hospital
# site code that appears in the synthetic registry+HES table.
#
# Columns the Stata step uses:
#   Trust_Name, Trust_Name_colour, Hospital_site_code, Bowel_ca_surgery,
#   Comprehensive_centre, Teaching_hospitals, Latest_Rating,
#   Staff_engagement, Moral, mean (= bed occupancy rate)
#
# Trust-level fields are held constant within a trust (the first 3 chars of the
# site code, matching the NHS nesting in the generator). Site-level fields vary
# by site. No real data is used.
# =============================================================================

library(tidyverse)
library(haven)
library(writexl)

base_dir <- "Data/synthetic/"
set.seed(20260601)

# site codes come straight from the synthetic registry+HES backbone
sites <- read_dta(paste0(base_dir, "colon_ncras_hes_synthetic.dta")) %>%
  distinct(diag_hosp) %>%
  filter(!is.na(diag_hosp)) %>%
  transmute(Hospital_site_code = diag_hosp,
            trust = substr(diag_hosp, 1, 3))   # site nests within trust prefix

# trust-level attributes, constant within a trust
trusts <- sites %>%
  distinct(trust) %>%
  mutate(
    Trust_Name = sprintf("Synthetic NHS Trust %03d", row_number()),
    # a minority of trusts get an excluded colour, to exercise the Stata drop
    Trust_Name_colour = sample(
      c("Blank","Yellow","Green","Grey","Light Red","Pink Red","Orange"),
      n(), replace = TRUE, prob = c(.45,.20,.15,.08,.05,.04,.03)),
    Staff_engagement = round(rnorm(n(), 6.90, 0.12), 4),
    Moral            = round(rnorm(n(), 5.95, 0.15), 4),
    mean             = round(pmin(pmax(rnorm(n(), 0.93, 0.025), 0.82), 0.99), 4)
  )

# site-level attributes joined to their trust
provider <- sites %>%
  left_join(trusts, by = "trust") %>%
  mutate(
    Bowel_ca_surgery     = sample(c(1, 0, NA), n(), replace = TRUE, prob = c(.55,.35,.10)),
    Comprehensive_centre = rbinom(n(), 1, 0.15),
    Teaching_hospitals   = rbinom(n(), 1, 0.25),
    Latest_Rating = sample(
      c("Outstanding","Good","Requires Improvement","Inadequate","Not rated"),
      n(), replace = TRUE, prob = c(.05,.45,.35,.08,.07)),
    mean = if_else(runif(n()) < 0.05, NA_real_, mean)   # a few missing, as in real
  ) %>%
  select(Trust_Name, Trust_Name_colour, Hospital_site_code,
         Bowel_ca_surgery, Comprehensive_centre, Teaching_hospitals,
         Latest_Rating, Staff_engagement, Moral, mean)

out_xlsx <- paste0(base_dir, "NHSHospitals_services_colon_synthetic.xlsx")
write_xlsx(provider, out_xlsx)

cat("Wrote", nrow(provider), "site rows to:\n  ", out_xlsx, "\n")
cat("Distinct trusts:", n_distinct(provider$Trust_Name), "\n")
cat("Colour-excluded sites:",
    sum(provider$Trust_Name_colour %in% c("Light Red","Pink Red","Orange")), "\n")
print(count(provider, Latest_Rating))
