# 03_travel_time_matrix_test_run.R
# computes car travel time matrix from lsoa centroids to hospital sites
# uses r5r with england osm network
# requires england-latest.osm.pbf in Data/travel_time/routing/

library(r5r)
library(dplyr)
library(purrr)
library(data.table)

dir.create("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/travel_time/routing/")
r5r_network <- build_network(data_path = "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/travel_time/routing/")

start <- Sys.time()

test_ttm <- travel_time_matrix(
  r5r_network,
  origins            = origins[1:10, ],
  destinations       = destinations[1:10, ],
  mode               = "CAR",
  departure_datetime = as.POSIXct("2024-06-11 10:00:00"),
  max_trip_duration  = 240L
)

cat("test time:", round(difftime(Sys.time(), start, units = "secs"), 1), "seconds\n")
cat("test rows:", nrow(test_ttm), "\n")
head(test_ttm)