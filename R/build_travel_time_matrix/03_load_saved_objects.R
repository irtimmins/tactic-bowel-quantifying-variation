# 03_load_saved_objects.R
# loads pre-processed objects from rds files
# run this instead of 01 and 02 in subsequent sessions

options(java.parameters = "-Xmx12G")

library(readxl)
library(data.table)
library(dplyr)
library(readr)
library(stringr)
library(sf)
library(purrr)
library(httr2)
library(r5r)

rds_path <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/rds"

hospital_df_postcodes <- readRDS(file.path(rds_path, "hospital_df_postcodes.rds"))
hospital_points       <- readRDS(file.path(rds_path, "hospital_points.rds"))
origins               <- readRDS(file.path(rds_path, "origins.rds"))
destinations          <- readRDS(file.path(rds_path, "destinations.rds"))
ods_lookup            <- readRDS(file.path(rds_path, "ods_lookup.rds"))
site_code_vec         <- readRDS(file.path(rds_path, "site_code_vec.rds"))
lsoa_centroids        <- readRDS(file.path(rds_path, "lsoa_centroids.rds"))

cat("loaded from", rds_path, "\n")
cat("origins:     ", nrow(origins), "\n")
cat("destinations:", nrow(destinations), "\n")
cat("site codes:  ", length(site_code_vec), "\n")