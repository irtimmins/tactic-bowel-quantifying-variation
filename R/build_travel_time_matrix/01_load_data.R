# 01_load_data.R
# loads raw data files needed for travel time matrix

options(java.parameters = "-Xmx12G")   # must be before r5r loads

library(readxl)
library(data.table)
library(dplyr)
library(readr)
library(stringr)
library(sf)
library(purrr)
library(httr2)
library(r5r)

base_path <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data"

# hospital curated dataset
hospital_df <- read_excel(file.path(base_path, "NHSHospitals_services_5.3.26_with_colours.xlsx"))

# site codes in scope, excluding pink/light red, dropping known duplicate
site_code_vec <- hospital_df %>%
  filter(trust_nacs_colour != "Pink Red", trust_nacs_colour != "Light Red") %>%
  filter(!is.na(Hospital_site_code)) %>%
  pull(Hospital_site_code) %>%
  unique() %>%
  .[. != "RN5T1"]   # duplicate of RN541 (royal hampshire winchester); ODS reused RN5T1 for basingstoke

cat("site codes in scope:", length(site_code_vec), "\n")

# lsoa population-weighted centroids (england and wales, 2011)
lsoa_centroids <- read_csv(
  file.path(base_path, "LSOA_Dec_2011_PWC_in_England_and_Wales_2022_1923591000694358693.csv")
) %>%
  select(lsoa11_code = LSOA11CD, easting = x, northing = y)

cat("lsoa centroids:", nrow(lsoa_centroids), "\n")

# onspd postcode to easting/northing lookup
onspd <- fread(file.path(base_path, "ONSPD_MAY_2025_UK.csv")) %>%
  select(postcode = pcds, easting = oseast1m, northing = osnrth1m) %>%
  filter(!is.na(easting))

cat("onspd postcodes:", nrow(onspd), "\n")

# ods site file
ods <- read_csv(file.path(base_path, "ets.csv"), col_names = FALSE) %>%
  select(site_code = X1, site_name = X2, postcode = X10,
         trust_code = X15, open_date = X11, close_date = X12) %>%
  mutate(
    open_date  = as.Date(as.character(open_date),  format = "%Y%m%d"),
    close_date = as.Date(as.character(close_date), format = "%Y%m%d")
  )

cat("ods sites loaded:", nrow(ods), "\n")