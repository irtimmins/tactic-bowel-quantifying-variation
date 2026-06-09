# 02_prepare_hospitals.R
# resolves hospital site codes to coordinates
# validates postcodes against ods api and applies corrections

library(dplyr)
library(stringr)
library(purrr)
library(httr2)
library(sf)

# helper functions
norm_pc <- function(x) toupper(gsub("\\s+", "", x))

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

to_lonlat <- function(df, id_col) {
  df %>%
    st_as_sf(coords = c("easting", "northing"), crs = 27700) %>%
    st_transform(4326) %>%
    mutate(lon = st_coordinates(.)[, 1],
           lat = st_coordinates(.)[, 2]) %>%
    st_drop_geometry() %>%
    transmute(id = .data[[id_col]], lon, lat)
}

# sites where the curated file has a typo or missing postcode
# RN5T1 already dropped from site_code_vec in 01_load_data.R
manual_postcodes <- tribble(
  ~site_code, ~site_name,                          ~postcode,
  "RTE26",    "Stroud General Hospital",            "GL5 2HY",   # na in curated file
  "RTE95",    "Tewkesbury Community Hospital",      "GL20 5QN",  # na in curated file
  "RBQ00",    "Liverpool Heart and Chest Hospital", "L14 3PE"    # not in ods
)

# build base postcode table from curated file
hospital_df_postcodes <- hospital_df %>%
  filter(trust_nacs_colour != "Pink Red", trust_nacs_colour != "Light Red") %>%
  filter(!is.na(Hospital_site_code)) %>%
  distinct(Hospital_site_code, .keep_all = TRUE) %>%
  select(site_code = Hospital_site_code,
         site_name = Hospital_Name,
         postcode  = Hospital_Post_Code) %>%
  mutate(postcode = case_when(
    site_code == "RAL16" ~ "WC1X 8DA",
    site_code == "RK5BC" ~ "NG17 4JL",
    TRUE ~ postcode
  )) %>%
  filter(site_code != "RN5T1") %>%      # duplicate of RN541, dropped from scope
  bind_rows(manual_postcodes) %>%
  filter(!is.na(postcode)) %>%
  distinct(site_code, .keep_all = TRUE)

cat("base postcodes built:", nrow(hospital_df_postcodes), "\n")

# query ods api for each site code to validate postcodes
fetch_ods_postcode <- function(code) {
  tryCatch({
    resp <- request("https://directory.spineservices.nhs.uk/ORD/2-0-0/organisations") %>%
      req_url_path_append(code) %>%
      req_perform()
    org <- resp_body_json(resp)$Organisation
    tibble(
      site_code  = code,
      ods_name   = org$Name,
      ods_status = org$Status,
      ods_pc     = org$GeoLoc$Location$PostCode %||% NA_character_
    )
  }, error = function(e) {
    tibble(site_code = code, ods_name = NA_character_,
           ods_status = "NOT_FOUND", ods_pc = NA_character_)
  })
}

cat("querying ods api for", length(site_code_vec), "site codes...\n")

ods_lookup <- map_dfr(site_code_vec, function(code) {
  Sys.sleep(0.15)
  fetch_ods_postcode(code)
}) %>%
  distinct(site_code, .keep_all = TRUE)

cat("ods api queries complete\n")

# apply ods postcodes, with named exceptions where curated is correct
hospital_df_postcodes <- hospital_df_postcodes %>%
  left_join(ods_lookup %>% select(site_code, ods_pc), by = "site_code") %>%
  mutate(postcode = case_when(
    site_code == "RBQ00" ~ "L14 3PE",   # not in ods, curated value confirmed
    site_code == "REN21" ~ "L9 7AL",    # clatterbridge liverpool -- ods correct
    site_code == "RRF01" ~ "WN7 1HS",   # leigh infirmary -- ods correct
    !is.na(ods_pc)       ~ ods_pc,      # prefer ods where available
    TRUE                 ~ postcode     # fall back to curated
  )) %>%
  select(-ods_pc) %>%
  distinct(site_code, .keep_all = TRUE)

# validation report
validation <- hospital_df_postcodes %>%
  left_join(ods_lookup %>% select(site_code, ods_name, ods_status, ods_pc),
            by = "site_code") %>%
  mutate(
    pc_curated = norm_pc(postcode),
    pc_ods     = norm_pc(ods_pc),
    match      = !is.na(pc_ods) & pc_curated == pc_ods
  ) %>%
  select(site_code, site_name, postcode, ods_pc, ods_name, ods_status, match)

cat("total sites:     ", nrow(validation), "\n")
cat("postcode matches:", sum(validation$match, na.rm = TRUE), "\n")
cat("mismatches:      ", sum(!validation$match, na.rm = TRUE), "\n")

remaining_mismatches <- validation %>% filter(!match)
if (nrow(remaining_mismatches) > 0) {
  cat("remaining mismatches:\n")
  print(remaining_mismatches)
}

# check all site_code_vec accounted for
missing <- setdiff(site_code_vec, hospital_df_postcodes$site_code)
if (length(missing) > 0) {
  cat("site codes with no postcode resolved:\n")
  print(missing)
} else {
  cat("all site codes resolved\n")
}

# join to onspd for coordinates
hospital_points <- hospital_df_postcodes %>%
  mutate(pc = norm_pc(postcode)) %>%
  left_join(onspd %>% mutate(pc = norm_pc(postcode)) %>%
              rename(onspd_postcode = postcode),
            by = "pc") %>%
  filter(!is.na(easting)) %>%
  distinct(site_code, .keep_all = TRUE)

cat("hospital points with coordinates:", nrow(hospital_points), "\n")

# any sites that failed the onspd join
no_coords <- setdiff(hospital_df_postcodes$site_code, hospital_points$site_code)
if (length(no_coords) > 0) {
  cat("sites with no coordinates:\n")
  print(hospital_df_postcodes %>% filter(site_code %in% no_coords))
}

# convert bng (27700) to lon/lat (4326) for r5r
origins      <- to_lonlat(lsoa_centroids, "lsoa11_code")
destinations <- to_lonlat(hospital_points, "site_code")

cat("origins:     ", nrow(origins), "\n")
cat("destinations:", nrow(destinations), "\n")