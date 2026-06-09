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

############################################################################


