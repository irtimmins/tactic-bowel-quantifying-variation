# 04_travel_time_matrix_full_run.R
# computes car travel time matrix from lsoa centroids to hospital sites
# uses r5r with england osm network
# requires england-latest.osm.pbf in Data/travel_time/routing/

library(r5r)
library(dplyr)
library(purrr)
library(data.table)

# recode lookup for downstream datasets -- anywhere RN5T1 appears recode to RN541
site_recode <- c("RN5T1" = "RN541")

# build network from osm pbf -- cached after first run as network.dat
dir.create("Data/travel_time/routing/")
r5r_network <- build_network(data_path = "Data/travel_time/routing/")

# split origins into chunks to keep java memory stable
chunks <- split(origins, (seq_len(nrow(origins)) - 1) %/% 3000)

cat("origins:     ", nrow(origins), "\n")
cat("destinations:", nrow(destinations), "\n")
cat("running matrix in", length(chunks), "chunks\n")

# compute matrix
ttm <- map_dfr(seq_along(chunks), function(i) {
  cat("chunk", i, "of", length(chunks), "\n")
  travel_time_matrix(
    r5r_network,
    origins            = chunks[[i]],
    destinations       = destinations,
    mode               = "CAR",
    departure_datetime = as.POSIXct("2024-10-15 08:30:00"),
    max_trip_duration  = 180L
  )
})

stop_r5(r5r_network)
rJava::.jgc(R.gc = TRUE)

cat("matrix rows:", nrow(ttm), "\n")

# nearest hospital per lsoa
nearest <- ttm %>%
  group_by(from_id) %>%
  slice_min(travel_time_p50, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(lsoa11_code     = from_id,
         nearest_site    = to_id,
         travel_time_min = travel_time_p50)

cat("lsoas with a nearest hospital:", nrow(nearest), "\n")

# lsoas with no hospital reachable within 180 min
missing_lsoas <- setdiff(origins$id, nearest$lsoa11_code)
if (length(missing_lsoas) > 0) {
  cat("lsoas with no reachable hospital:", length(missing_lsoas), "\n")
  print(head(missing_lsoas))
} else {
  cat("all lsoas have at least one reachable hospital\n")
}

# save outputs
fwrite(ttm,     "Data/travel_time/ttm_lsoa_hospital.csv")
fwrite(nearest, "Data/travel_time/nearest_hospital_per_lsoa.csv")

cat("done\n")