library(haven)

lu_matrix <- read_dta("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/travel_time/bowel_pairwise_distance_matrix_lu.dta")

# sites in colleagues matrix
lu_sites <- lu_matrix %>%
  pull(sitecode) %>%
  unique() %>%
  sort()

cat("sites in colleagues matrix:", length(lu_sites), "\n")
cat("sites in your matrix:      ", length(site_code_vec), "\n")

# what's in yours but not theirs
only_in_yours <- setdiff(site_code_vec, lu_sites)
cat("\nin yours but not colleagues (", length(only_in_yours), "):\n")
print(only_in_yours)

# what's in theirs but not yours
only_in_theirs <- setdiff(lu_sites, site_code_vec)
cat("\nin colleagues but not yours (", length(only_in_theirs), "):\n")
#print(only_in_theirs)
# check what these 13 are
ods %>%
  filter(site_code %in% only_in_theirs) %>%
  distinct(site_code, .keep_all = TRUE) %>%
  select(site_code, site_name, postcode, open_date, close_date) %>%
  arrange(site_code)

# overlap
both <- intersect(site_code_vec, lu_sites)
cat("\nin both:", length(both), "\n")

# lsoa coverage check
lu_lsoas <- lu_matrix %>% pull(lsoa11_code) %>% unique()
cat("\nlsoas in colleagues matrix:", length(lu_lsoas), "\n")
cat("lsoas in your origins:     ", nrow(origins), "\n")

# do any of the 13 postcodes match postcodes in your hospital_df?
theirs_postcodes <- ods %>%
  filter(site_code %in% only_in_theirs) %>%
  distinct(site_code, site_name, postcode)


hospital_df_postcodes %>%
  mutate(pc = norm_pc(postcode)) %>%
  inner_join(theirs_postcodes %>%
               mutate(pc = norm_pc(postcode)) %>%
               select(site_code_theirs = site_code,
                      site_name_theirs = site_name,
                      pc),
             by = "pc") %>%
  select(your_code = site_code, your_name = site_name,
         their_code = site_code_theirs, their_name = site_name_theirs, postcode)

# what are the 79 sites in yours but not colleagues
yours_only_detail <- map_dfr(only_in_yours, function(code) {
  Sys.sleep(0.15)
  fetch_ods_postcode(code)
}) %>%
  select(site_code, ods_name, ods_status, ods_pc)

print(yours_only_detail, n = 79)

hospital_df %>%
  filter(Hospital_site_code %in% only_in_yours) %>%
  distinct(Hospital_site_code, .keep_all = TRUE) %>%
  select(site_code = Hospital_site_code,
         site_name = Hospital_Name,
         trust_name = Trust_Name,
         trust_colour = trust_nacs_colour) %>%
  left_join(yours_only_detail, by = "site_code") %>%
  arrange(trust_colour, site_code) %>%
  print(n = 79)


# investigate the 1869 lsoas missing from your matrix
cat("lsoas in lu but not yours:\n")
missing_lsoas_detail <- lu_nearest %>%
  filter(lsoa11_code %in% only_lu) %>%
  summarise(
    n             = n(),
    min_lu_time   = min(travel_time_lu),
    median_lu_time = median(travel_time_lu),
    max_lu_time   = max(travel_time_lu)
  )
print(missing_lsoas_detail)

# are they wales lsoas?
lu_nearest %>%
  filter(lsoa11_code %in% only_lu) %>%
  mutate(country = substr(lsoa11_code, 1, 1)) %>%
  count(country)


# investigate E01011572 -- same hospital, 10x time difference
# check what's in the full matrix for this lsoa
ttm_full <- fread(file.path(output_path, "ttm_lsoa_hospital.csv"))

ttm_full %>%
  filter(from_id == "E01011572") %>%
  arrange(travel_time_p50) %>%
  head(10)

# compare with lu
ttm_lu %>%
  filter(lsoa11_code == "E01011572") %>%
  arrange(total_drive_time) %>%
  head(10)



# look at the rnlbx pattern more carefully
# these are correct in your matrix -- lu simply didn't have rnlbx
rnlbx_check <- comparison %>%
  filter(nearest_site == "RNLBX" | nearest_site_lu == "RNLAY") %>%
  select(lsoa11_code, nearest_site, nearest_site_lu,
         travel_time_min, travel_time_lu, time_diff)

cat("rnlbx/rnlay discrepancies:", nrow(rnlbx_check), "\n")
print(head(rnlbx_check, 10))


# rerun comparison excluding known hospital list differences
# recode lu matrix sites to your codes where we know they match
lu_recoded <- lu_nearest %>%
  mutate(nearest_site_lu = recode(nearest_site_lu,
    "RNLAY"  = "RNLBX",   # same physical site
    "RC110"  = "RC979",
    "RGQ02"  = "RDE03",
    "RQ617"  = "REMRQ",
    "RDDH0"  = "RAJ12",
    "RQ8L0"  = "RAJ32",
    "R0B0Q"  = "RX454",
    "RE9GA"  = "RX454",
    "R0B01"  = "RX4J5",
    "RLNGL"  = "RX4J5",
    "RR105"  = "RRK98",
    "RR101"  = "RRK97",
    "RA301"  = "RA7C2",
    "RJF02"  = "RTG02"
  ))

comparison_recoded <- ttm %>%
  inner_join(lu_recoded, by = "lsoa11_code") %>%
  mutate(time_diff = travel_time_min - travel_time_lu,
         same      = nearest_site == nearest_site_lu)

cat("\nafter recoding same-site codes:\n")
comparison_recoded %>%
  count(same) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

cat("\ncorrelation after recoding:", 
    round(cor(comparison_recoded$travel_time_min,
              comparison_recoded$travel_time_lu), 4), "\n")


############################################################################


# find the one missing english lsoa
missing_english <- lu_nearest %>%
  filter(lsoa11_code %in% only_lu,
         substr(lsoa11_code, 1, 1) == "E") %>%
  pull(lsoa11_code)
cat("missing english lsoa:", missing_english, "\n")

# check if it's in origins
missing_english %in% origins$id

# look at disagreements more carefully
# split by whether same or different nearest hospital
cat("\n--- time diff where same nearest hospital ---\n")
comparison_recoded %>%
  filter(same) %>%
  summarise(
    n          = n(),
    mean_diff  = round(mean(time_diff), 2),
    median_diff = round(median(time_diff), 2),
    sd_diff    = round(sd(time_diff), 2),
    pct_within_5  = round(100 * mean(abs(time_diff) < 5), 1),
    pct_within_10 = round(100 * mean(abs(time_diff) < 10), 1),
    pct_within_15 = round(100 * mean(abs(time_diff) < 15), 1)
  ) %>%
  print()

cat("\n--- time diff where different nearest hospital ---\n")
comparison_recoded %>%
  filter(!same) %>%
  summarise(
    n          = n(),
    mean_diff  = round(mean(time_diff), 2),
    median_diff = round(median(time_diff), 2),
    sd_diff    = round(sd(time_diff), 2)
  ) %>%
  print()

# where same hospital -- are differences systematic?
# r5r returns integer minutes, lu is continuous
# so expected difference is 0-1 min from rounding alone
cat("\n--- distribution of diffs where same hospital ---\n")
comparison_recoded %>%
  filter(same) %>%
  mutate(diff_band = cut(time_diff,
                         breaks = c(-Inf, -10, -5, 0, 5, 10, 20, Inf),
                         labels = c("<-10", "-10 to -5", "-5 to 0",
                                    "0 to 5", "5 to 10", "10 to 20", ">20"))) %>%
  count(diff_band) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

# where different hospital -- is yours systematically shorter?
cat("\n--- direction of difference where different hospital ---\n")
comparison_recoded %>%
  filter(!same) %>%
  mutate(direction = case_when(
    time_diff < -5  ~ "yours much shorter",
    time_diff < 0   ~ "yours slightly shorter",
    time_diff < 5   ~ "similar",
    time_diff < 10  ~ "yours slightly longer",
    TRUE            ~ "yours much longer"
  )) %>%
  count(direction) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n)) %>%
  print()

# check E01011572 -- is it an island or remote lsoa?
cat("\nE01011572 in origins:\n")
origins %>% filter(id == "E01011572")


##############################################################



origins %>% filter(id == "E01028883")



# check E01028883
origins %>% filter(id == "E01028883")

# check what lu has for E01028883
ttm_lu %>%
  filter(lsoa11_code == "E01028883") %>%
  arrange(total_drive_time) %>%
  head(5)

# investigate E01011572 more -- what's nearby in the road network?
# check neighbouring lsoas for comparison
lsoa_centroids %>%
  filter(lsoa11_code %in% c("E01011572", "E01011571", "E01011573",
                             "E01011574", "E01011575")) %>%
  print()

# check if E01011572 times are all inflated or just some
ttm_full %>%
  filter(from_id == "E01011572") %>%
  arrange(travel_time_p50) %>%
  print()

# check neighbouring lsoas travel times to RAE01 for context
ttm_full %>%
  filter(from_id %in% c("E01011571", "E01011573", "E01011574"),
         to_id == "RAE01") %>%
  arrange(from_id)



############################################################
######################

# identify all potential snapping failures
# flag lsoas where travel time to nearest hospital is >3x the median
# of their geographic neighbours

# first get a sense of scale -- how many lsoas have suspiciously long times
ttm %>%
  mutate(suspicious = travel_time_min > 120) %>%
  count(suspicious) %>%
  mutate(pct = round(100 * n / sum(n), 1))



# compare each lsoa's nearest time against lu's nearest time
# large positive differences (yours >> lu) suggest snapping issues
snapping_suspects <- comparison_recoded %>%
  filter(same,
         time_diff > 30) %>%
  select(lsoa11_code, nearest_site, travel_time_min,
         travel_time_lu, time_diff) %>%
  arrange(desc(time_diff))

cat("potential snapping failures (same hospital, >30 min difference):", 
    nrow(snapping_suspects), "\n")
print(snapping_suspects)

# for the methodology section -- overall quality summary
cat("=== matrix quality summary ===\n")
cat("total lsoas in matrix:          ", nrow(ttm), "\n")
cat("lsoas with no result:           ", 1, "(E01028883)\n")
cat("potential snapping failures:    ", nrow(snapping_suspects), "\n")
cat("agreement with lu (same hosp):  ",
    round(100 * mean(comparison_recoded$same), 1), "%\n")
cat("within 10 min where same hosp:  97.1%\n")
cat("systematic offset vs lu:        +2.3 mins (routing engine difference)\n")
cat("welsh lsoas excluded:           1868 (england-only analysis)\n")

###########################################################################
####################################################


# check E01009543 neighbours
ttm_full %>%
  filter(from_id %in% c("E01009543", "E01009542", "E01009544", "E01009545"),
         to_id == "RKB01") %>%
  arrange(from_id)

origins %>% filter(id == "E01009543")

# fix snapping failures by substituting neighbour median
# E01011572 -- neighbours E01011571, E01011573, E01011574
# E01009543 -- check neighbours first then fix

fix_lsoas <- c("E01011572", "E01009543", "E01028883")

# get all results for neighbours of E01011572
neighbour_times_11572 <- ttm_full %>%
  filter(from_id %in% c("E01011571", "E01011573", "E01011574")) %>%
  group_by(to_id) %>%
  summarise(travel_time_p50 = round(median(travel_time_p50)),
            .groups = "drop") %>%
  mutate(from_id = "E01011572")

# get all results for neighbours of E01009543
neighbour_times_09543 <- ttm_full %>%
  filter(from_id %in% c("E01009542", "E01009544", "E01009545")) %>%
  group_by(to_id) %>%
  summarise(travel_time_p50 = round(median(travel_time_p50)),
            .groups = "drop") %>%
  mutate(from_id = "E01009543")

# get all results for neighbours of E01028883 -- use lu times as proxy
neighbour_times_28883 <- ttm_lu %>%
  filter(lsoa11_code == "E01028883") %>%
  transmute(from_id        = "E01028883",
            to_id          = sitecode,
            travel_time_p50 = as.integer(round(total_drive_time)))

# remove failed rows and add fixed rows
ttm_fixed <- ttm_full %>%
  filter(!from_id %in% fix_lsoas) %>%
  bind_rows(neighbour_times_11572,
            neighbour_times_09543,
            neighbour_times_28883)

cat("rows before fix:", nrow(ttm_full), "\n")
cat("rows after fix: ", nrow(ttm_fixed), "\n")
cat("unique origins: ", n_distinct(ttm_fixed$from_id), "\n")

# rebuild nearest hospital
nearest_fixed <- ttm_fixed %>%
  group_by(from_id) %>%
  slice_min(travel_time_p50, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(lsoa11_code     = from_id,
         nearest_site    = to_id,
         travel_time_min = travel_time_p50)

cat("lsoas with nearest hospital:", nrow(nearest_fixed), "\n")

# save fixed versions
fwrite(ttm_fixed,     file.path(output_path, "ttm_lsoa_hospital_fixed.csv"))
fwrite(nearest_fixed, file.path(output_path, "nearest_hospital_per_lsoa_fixed.csv"))

cat("saved fixed matrix\n")


##########################################
#############

cat("=== final matrix summary ===\n")
cat("total lsoa-hospital pairs:   ", nrow(ttm_fixed), "\n")
cat("unique lsoas:                ", n_distinct(ttm_fixed$from_id), "\n")
cat("unique hospitals:            ", n_distinct(ttm_fixed$to_id), "\n")
cat("travel time range (mins):    ", 
    min(ttm_fixed$travel_time_p50), "to",
    max(ttm_fixed$travel_time_p50), "\n")

cat("\n--- nearest hospital distribution ---\n")
nearest_fixed %>%
  mutate(time_band = cut(travel_time_min,
                         breaks = c(0, 15, 30, 45, 60, 90, 120, Inf),
                         labels = c("<15", "15-30", "30-45",
                                    "45-60", "60-90", "90-120", ">120"))) %>%
  count(time_band) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

cat("\n--- corrections applied ---\n")
cat("E01011572: snapping failure, imputed from 3 neighbours\n")
cat("E01009543: snapping failure, imputed from 2 neighbours\n")
cat("E01028883: no route found, imputed from lu reference matrix\n")

################################################
###########

# hospitals in final matrix
hospitals_in_matrix <- ttm_fixed %>%
  pull(to_id) %>%
  unique() %>%
  sort()

cat("hospitals in final matrix:", length(hospitals_in_matrix), "\n")
cat("hospitals in scope:       ", length(site_code_vec), "\n")

# missing from matrix
missing_hospitals <- setdiff(site_code_vec, hospitals_in_matrix)
cat("missing from matrix:      ", length(missing_hospitals), "\n")
print(missing_hospitals)

# in matrix but not in scope
extra_hospitals <- setdiff(hospitals_in_matrix, site_code_vec)
cat("in matrix but not in scope:", length(extra_hospitals), "\n")
print(extra_hospitals)


#########################################################


# remove legacy codes that crept in via lu imputation
ttm_fixed <- ttm_fixed %>%
  filter(!to_id %in% extra_hospitals)

# rebuild nearest
nearest_fixed <- ttm_fixed %>%
  group_by(from_id) %>%
  slice_min(travel_time_p50, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(lsoa11_code     = from_id,
         nearest_site    = to_id,
         travel_time_min = travel_time_p50)

# confirm
cat("hospitals in fixed matrix:", n_distinct(ttm_fixed$to_id), "\n")
cat("unique lsoas:             ", n_distinct(ttm_fixed$from_id), "\n")
cat("total rows:               ", nrow(ttm_fixed), "\n")
cat("missing hospitals:        ", 
    length(setdiff(site_code_vec, unique(ttm_fixed$to_id))), "\n")

# resave
fwrite(ttm_fixed,     file.path(output_path, "ttm_lsoa_hospital_fixed.csv"))
fwrite(nearest_fixed, file.path(output_path, "nearest_hospital_per_lsoa_fixed.csv"))

cat("saved\n")

##############################################################



# hospitals in scope that lu's matrix didn't have
lu_sites <- ttm_lu %>% pull(sitecode) %>% unique()

in_yours_not_lu <- setdiff(site_code_vec, lu_sites)
in_both          <- intersect(site_code_vec, lu_sites)
in_lu_not_yours  <- setdiff(lu_sites, site_code_vec)

cat("hospitals in your matrix:        ", length(site_code_vec), "\n")
cat("hospitals in lu matrix:          ", length(lu_sites), "\n")
cat("in both:                         ", length(in_both), "\n")
cat("in yours but not lu:             ", length(in_yours_not_lu), "\n")
cat("in lu but not yours (legacy):    ", length(in_lu_not_yours), "\n")

cat("\nhospitals lu missed:\n")
hospital_df_postcodes %>%
  filter(site_code %in% in_yours_not_lu) %>%
  select(site_code, site_name, postcode) %>%
  arrange(site_name) %>%
  print(n = 100)

###############################################################

# hospital curated dataset
hospital_df <- read_excel(file.path(base_path, "NHSHospitals_services_5.3.26_with_colours.xlsx"))


# any pattern in what lu missed -- trust level?
hospital_df %>%
  filter(Hospital_site_code %in% in_yours_not_lu) %>%
  distinct(Hospital_site_code, .keep_all = TRUE) %>%
  count(trust_nacs_colour) %>%
  arrange(desc(n))


# any pattern by service type?
hospital_df %>%
  filter(Hospital_site_code %in% in_yours_not_lu) %>%
  distinct(Hospital_site_code, .keep_all = TRUE) %>%
  summarise(
    bowel_surgery    = sum(Bowel_ca_surgery    == "Y", na.rm = TRUE),
    lung_surgery     = sum(Lung_Ca_surgery     == "Y", na.rm = TRUE),
    radiotherapy     = sum(Radiotherapy        == "Y", na.rm = TRUE),
    chemotherapy     = sum(Chemo               == "Y", na.rm = TRUE),
    hepatobiliary    = sum(Hepatobiliary_surgery == "Y", na.rm = TRUE),
    comprehensive    = sum(Comprehensive_centre == "Y", na.rm = TRUE)
  ) %>%
  print()



# note two hospitals appear twice -- broadgreen and queen marys sidcup
# check why
hospital_df_postcodes %>%
  filter(site_code %in% c("REMAH", "RQ601", "RN7QM", "RJ230"))



####################################################

# confirm they have identical travel times
ttm_fixed %>%
  filter(to_id %in% c("RN7QM", "RJ230")) %>%
  group_by(to_id) %>%
  summarise(mean_time = mean(travel_time_p50),
            n = n()) %>%
  print()

ttm_fixed %>%
  filter(to_id %in% c("REMAH", "RQ601")) %>%
  group_by(to_id) %>%
  summarise(mean_time = mean(travel_time_p50),
            n = n()) %>%
  print()


###
# check what these lsoas and hospitals are
hospital_df_postcodes %>%
  filter(site_code %in% c("RJ611", "RYJ01")) %>%
  select(site_code, site_name, postcode)

origins %>%
  filter(id %in% c("E01001178", "E01033594"))

# check neighbouring lsoas travel time to same hospital for context
ttm_fixed %>%
  filter(to_id == "RJ611",
         from_id %in% c("E01001177", "E01001178", "E01001179", "E01001180")) %>%
  arrange(travel_time_p50)

ttm_fixed %>%
  filter(to_id == "RYJ01",
         from_id %in% c("E01033593", "E01033594", "E01033595", "E01033596")) %>%
  arrange(travel_time_p50)


########################################

# quick export to stata
library(haven)
ttm_final <- ttm_fixed %>%
            rename(lsoa11_code     = from_id,
                   site_code       = to_id,
                   travel_time_min = travel_time_p50)
write_dta(ttm_final,
          file.path(output_path, "travel_time_matrix_tactic_260611.dta"))

write_dta(nearest_fixed,
          file.path(output_path, "nearest_hospital_per_lsoa_fixed.dta"))

cat("exported to stata format\n")
















