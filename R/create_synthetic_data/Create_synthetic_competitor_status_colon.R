# =============================================================================
# Synthetic competitor status (net patient flow) for the colon provider analysis
# -----------------------------------------------------------------------------
# Classifies each hospital site as a net importer ("Winner"), net exporter
# ("Loser") or no significant difference, and writes the file the Stata
# provider step optionally merges:
#
#   colon_competitor_status_synthetic.dta -> diag_hosp, competitor_status
#       competitor_status: 1 = Winner, 2 = Loser, 3 = Insignificant diff.
#
# The colon registry+HES table does not carry an LSOA of residence, so the
# patient geography is synthesised here: each site gets a 2-D location and an
# "attractiveness" (pull); each patient is given a home district at the lowest
# pull site among co-attenders, and a patient diagnosed at a higher-pull site
# counts as an arriver there and a leaver at the home site. High-pull sites
# therefore import and low-pull sites export, with leavers and arrivers
# balancing overall. No real data is used.
#
# It also writes the distance matrix and a valid-sites lookup, which are not
# needed by the current Stata step but are kept for parity with the real
# win/loss workflow.
# =============================================================================

library(tidyverse)
library(haven)

base_dir <- "Data/synthetic/"
set.seed(20260601)

# ---- realism knobs ----------------------------------------------------------
SPACE       <- 100     # size of the synthetic map (arbitrary units)
SCALE_MIN   <- 2.0     # drive-time minutes per map unit
HOME_SD     <- 1.0     # scatter of a district around its local site
NOISE_SD    <- 0.3     # noise on each drive time (minutes)
PULL_SPREAD <- 0.8     # >1 sharpens the winner/loser contrast, <1 softens it
N_DISTRICT  <- 250L    # number of synthetic residence districts
SIG_CUT     <- 1.96    # net-flow z-score for a significant winner/loser

# ---- patients and their diagnosing site -------------------------------------
cohort <- read_dta(paste0(base_dir, "colon_ncras_hes_synthetic.dta")) %>%
  select(pseudo_patientid, diag_hosp) %>%
  filter(!is.na(diag_hosp))

# ---- site geography: coordinates and attractiveness -------------------------
sites <- tibble(sitecode = sort(unique(cohort$diag_hosp))) %>%
  mutate(sx   = runif(n(), 0, SPACE),
         sy   = runif(n(), 0, SPACE),
         pull = (rank(runif(n())) / (n() + 1))^(1 / PULL_SPREAD))
cat("Synthetic hospital sites:", nrow(sites), "\n")

# ---- give each patient a synthetic home district ----------------------------
# districts sit near a randomly chosen "local" site; patients are assigned to a
# district, so co-resident patients can attend different (higher-pull) sites
districts <- tibble(
  lsoa11_code = sprintf("E0%06d", seq_len(N_DISTRICT)),   # synthetic LSOA codes
  local_site = sample(sites$sitecode, N_DISTRICT, replace = TRUE,
                      prob = 1 - sites$pull[match(sites$sitecode, sites$sitecode)] + 0.1)
) %>%
  left_join(sites, by = c("local_site" = "sitecode")) %>%
  transmute(lsoa11_code, local_site,
            hx = sx + rnorm(n(), 0, HOME_SD),
            hy = sy + rnorm(n(), 0, HOME_SD))

pat <- cohort %>%
  mutate(lsoa11_code = sample(districts$lsoa11_code, n(), replace = TRUE)) %>%
  left_join(districts %>% select(lsoa11_code, hx, hy), by = "lsoa11_code")

# ---- nearest site to each lsoa (the "expected" local provider) --------------
dist_long <- crossing(districts %>% select(lsoa11_code, hx, hy),
                      sites %>% select(sitecode, sx, sy)) %>%
  mutate(total_drive_time = round(
    sqrt((hx - sx)^2 + (hy - sy)^2) * SCALE_MIN + 1 + abs(rnorm(n(), 0, NOISE_SD)), 2))

nearest <- dist_long %>%
  group_by(lsoa11_code) %>%
  slice_min(total_drive_time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(lsoa11_code, nearest_site = sitecode)

# ---- net patient flow per site ----------------------------------------------
flows <- pat %>%
  left_join(nearest, by = "lsoa11_code") %>%
  mutate(core = diag_hosp == nearest_site)

leavers  <- flows %>% filter(!core) %>% count(site = nearest_site, name = "n_leavers")
arrivers <- flows %>% filter(!core) %>% count(site = diag_hosp,    name = "n_arrivers")

site_flow <- sites %>%
  transmute(site = sitecode, pull) %>%
  left_join(leavers,  by = "site") %>%
  left_join(arrivers, by = "site") %>%
  mutate(across(c(n_leavers, n_arrivers), ~replace_na(.x, 0L)),
         n_net = n_arrivers - n_leavers,
         # a simple z-score on the net flow, treating arrivers+leavers as the base
         base  = pmax(n_arrivers + n_leavers, 1),
         z     = n_net / sqrt(base),
         competitor_status = case_when(
           z >=  SIG_CUT ~ 1L,   # Winner  (net importer)
           z <= -SIG_CUT ~ 2L,   # Loser   (net exporter)
           TRUE          ~ 3L))  # Insignificant difference

# ---- write the file the Stata step merges -----------------------------------
comp <- site_flow %>% transmute(diag_hosp = site, competitor_status)
write_dta(comp, paste0(base_dir, "colon_competitor_status_synthetic.dta"))

# patient -> lsoa11_code lookup: used by the Stata step on synthetic data
# (on real data lsoa11_code is already in the cohort)
pat_lsoa <- pat %>% select(pseudo_patientid, lsoa11_code)
write_dta(pat_lsoa, paste0(base_dir, "colon_patient_lsoa_synthetic.dta"))

# ---- also write distance matrix + valid sites, for parity with the real flow -
dist_out <- dist_long %>% transmute(lsoa11_code, sitecode, total_drive_time)
write_dta(dist_out,  paste0(base_dir, "colon_pairwise_distance_matrix_synthetic.dta"))
valid <- sites %>% transmute(diag_hosp = sitecode, valid = 1L)
write_dta(valid, paste0(base_dir, "colon_valid_sites_synthetic.dta"))

# ---- quick QC ---------------------------------------------------------------
cat("Total leavers:", sum(site_flow$n_leavers),
    "| total arrivers:", sum(site_flow$n_arrivers), "(should match)\n")
cat("Patients in a flow:", round(100 * mean(!flows$core), 1), "%\n")
print(count(comp, competitor_status))
cat("Correlation(pull, net flow):", round(cor(site_flow$pull, site_flow$n_net), 2), "\n")
