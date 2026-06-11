# benchmark_r5r.R
# systematic timing tests to estimate full run time
# tests different chunk sizes and max_trip_durations

library(rslurm)

source("R/build_travel_time_matrix/03_load_saved_objects.R")
output_path  <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/output"
routing_base <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/routing"

rds_path <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/rds"

if (!exists("origins")) {
  cat("loading saved objects\n")
  origins      <- readRDS(file.path(rds_path, "origins.rds"))
  destinations <- readRDS(file.path(rds_path, "destinations.rds"))
  cat("origins:     ", nrow(origins), "\n")
  cat("destinations:", nrow(destinations), "\n")
}

if (!exists("origin_chunks")) {
  n_chunks      <- 500
  chunk_id      <- (seq_len(nrow(origins)) - 1) %/% ceiling(nrow(origins) / n_chunks)
  origin_chunks <- split(origins, chunk_id)
  cat("chunks built:", length(origin_chunks), "\n")
}

benchmark_r5r <- function(chunk_index, n_origins, max_duration) {

  options(java.parameters = "-Xmx24G")
  library(r5r)
  library(data.table)

  routing_base <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/routing"
  output_path  <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/output"

  # isolated working directory
  job_dir <- file.path(tempdir(), paste0("bench_", chunk_index))
  dir.create(job_dir, showWarnings = FALSE)
  old_wd <- setwd(job_dir)
  on.exit(setwd(old_wd))

  # copy network locally
  local_routing <- file.path(job_dir, "routing")
  dir.create(local_routing, showWarnings = FALSE)

  t_copy_start <- Sys.time()
  file.copy(file.path(routing_base, "england-260608.osm.pbf"),
            file.path(local_routing, "england-260608.osm.pbf"))
  file.copy(file.path(routing_base, "network.dat"),
            file.path(local_routing, "network.dat"))
  t_copy <- as.numeric(difftime(Sys.time(), t_copy_start, units = "secs"))

  # build network
  t_build_start <- Sys.time()
  r5r_network <- build_network(data_path = local_routing, verbose = FALSE)
  t_build <- as.numeric(difftime(Sys.time(), t_build_start, units = "secs"))

  # use real origins and destinations
  origins_chunk <- origin_chunks[[chunk_index]][1:n_origins, ]

  # routing
  t_route_start <- Sys.time()
  ttm <- travel_time_matrix(
    r5r_network,
    origins            = origins_chunk,
    destinations       = destinations,
    mode               = "CAR",
    departure_datetime = as.POSIXct("2024-06-11 10:00:00"),
    max_trip_duration  = as.integer(max_duration)
  )
  t_route <- as.numeric(difftime(Sys.time(), t_route_start, units = "secs"))

  stop_r5(r5r_network)
  rJava::.jgc(R.gc = TRUE)
  unlink(job_dir, recursive = TRUE)

  # results
  result <- data.frame(
    chunk_index   = chunk_index,
    n_origins     = n_origins,
    max_duration  = max_duration,
    n_destinations = nrow(destinations),
    t_copy_secs   = round(t_copy,   1),
    t_build_secs  = round(t_build,  1),
    t_route_secs  = round(t_route,  1),
    t_total_secs  = round(t_copy + t_build + t_route, 1),
    n_pairs       = nrow(ttm),
    node          = Sys.info()["nodename"]
  )

  cat("chunk:", chunk_index,
      "| origins:", n_origins,
      "| max_dur:", max_duration,
      "| copy:", round(t_copy, 1), "s",
      "| build:", round(t_build, 1), "s",
      "| route:", round(t_route, 1), "s",
      "| total:", round(t_copy + t_build + t_route, 1), "s\n")

  # save individual result
  fwrite(result,
         file.path(output_path,
                   paste0("bench_", chunk_index, "_n", n_origins,
                          "_d", max_duration, ".csv")))

  return(result)
}

# benchmark parameters
# vary origins: 10, 50, 100, 200, 500
# vary max_duration: 120, 240, 360, 480, 560
# use different chunk_index values to get different geographic spread of origins
pars_bench <- expand.grid(
  chunk_index  = 1:5,           # 3 geographic areas
  n_origins    = c(10, 50, 100, 200),
  max_duration = c(120, 240, 360, 600)
) %>%
  arrange(chunk_index, n_origins, max_duration)

cat("total benchmark jobs:", nrow(pars_bench), "\n")
print(pars_bench)

sjob_bench <- slurm_apply(
  benchmark_r5r,
  pars_bench,
  jobname        = "r5r_benchv2",
  nodes          = nrow(pars_bench),
  cpus_per_node  = 1,
  submit         = FALSE,
  global_objects = c("origin_chunks", "destinations"),
  pkgs           = c("r5r", "data.table"),
  slurm_options  = list(
    time          = "06:00:00",
    partition     = "normal",
    "mem-per-cpu" = "32G"
  )
)

# patch
log_bench <- list.files(
  "/home/lshit9/TACTIC/build_distance_matrix",
  pattern    = "_rslurm_r5r_bench",
  full.names = TRUE
)

run_script <- readLines(file.path(log_bench, "slurm_run.R"))
writeLines(c('options(java.parameters = "-Xmx24G")', run_script),
           file.path(log_bench, "slurm_run.R"))

submit_script <- readLines(file.path(log_bench, "submit.sh"))
submit_script <- gsub("--vanilla", "", submit_script)
writeLines(submit_script, file.path(log_bench, "submit.sh"))

cat("--- first 3 lines slurm_run.R ---\n")
cat(head(readLines(file.path(log_bench, "slurm_run.R")), 3), sep = "\n")
cat("\n--- submit.sh ---\n")
cat(readLines(file.path(log_bench, "submit.sh")), sep = "\n")

# submit
system(paste("cd", log_bench, "&& sbatch submit.sh"))
