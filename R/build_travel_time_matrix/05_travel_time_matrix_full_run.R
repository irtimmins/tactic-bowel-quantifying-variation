# 05 
# full travel time matrix -- conservative settings
# 100 origins per chunk, 600 min max duration, 4h time limit
# requires 03_load_saved_objects.R to have been sourced first

source("/home/lshit9/TACTIC/build_distance_matrix/R/build_travel_time_matrix/03_load_saved_objects.R")


library(rslurm)
library(dplyr)
library(data.table)

base_path    <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data"
routing_base <- file.path(base_path, "routing")
output_path  <- file.path(base_path, "output")
rds_path     <- file.path(base_path, "rds")

dir.create(output_path, showWarnings = FALSE)

# load saved objects if not in session
if (!exists("origins")) {
  cat("loading saved objects\n")
  origins      <- readRDS(file.path(rds_path, "origins.rds"))
  destinations <- readRDS(file.path(rds_path, "destinations.rds"))
  cat("origins:     ", nrow(origins), "\n")
  cat("destinations:", nrow(destinations), "\n")
}

# build chunks
n_origins_per_chunk <- 100
chunk_id            <- (seq_len(nrow(origins)) - 1) %/% n_origins_per_chunk
origin_chunks       <- split(origins, chunk_id)
n_chunks            <- length(origin_chunks)

cat("origins per chunk:", n_origins_per_chunk, "\n")
cat("total chunks:     ", n_chunks, "\n")

# function to run on each node
run_ttm_chunk <- function(chunk_index) {

  options(java.parameters = "-Xmx48G")
  library(r5r)
  library(data.table)

  routing_base <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/routing"
  output_path  <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/output"

  # isolated working directory per job
  job_dir <- file.path(tempdir(), paste0("job_", chunk_index))
  dir.create(job_dir, showWarnings = FALSE)
  old_wd <- setwd(job_dir)
  on.exit({
    setwd(old_wd)
    unlink(job_dir, recursive = TRUE)
  })

  # copy network to local temp
  local_routing <- file.path(job_dir, "routing")
  dir.create(local_routing, showWarnings = FALSE)

  t_start <- Sys.time()

  file.copy(file.path(routing_base, "england-260608.osm.pbf"),
            file.path(local_routing, "england-260608.osm.pbf"))
  file.copy(file.path(routing_base, "network.dat"),
            file.path(local_routing, "network.dat"))

  cat(format(Sys.time(), "%H:%M:%S"), "chunk", chunk_index,
      "-- network copied\n")

  r5r_network <- build_network(data_path = local_routing, verbose = FALSE)

  cat(format(Sys.time(), "%H:%M:%S"), "chunk", chunk_index,
      "-- network built, routing", nrow(origin_chunks[[chunk_index]]),
      "origins x", nrow(destinations), "destinations\n")

  ttm_chunk <- travel_time_matrix(
    r5r_network,
    origins            = origin_chunks[[chunk_index]],
    destinations       = destinations,
    mode               = "CAR",
    departure_datetime = as.POSIXct("2024-06-11 10:00:00"),
    max_trip_duration  = 600L
  )

  stop_r5(r5r_network)
  rJava::.jgc(R.gc = TRUE)

  t_elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1)

  cat(format(Sys.time(), "%H:%M:%S"), "chunk", chunk_index,
      "-- done:", nrow(ttm_chunk), "rows,", t_elapsed, "mins\n")

  fwrite(ttm_chunk,
         file.path(output_path, paste0("ttm_chunk_", chunk_index, ".csv")))

  return(data.frame(
    chunk_index  = chunk_index,
    n_origins    = nrow(origin_chunks[[chunk_index]]),
    n_rows       = nrow(ttm_chunk),
    elapsed_min  = t_elapsed,
    node         = Sys.info()["nodename"]
  ))
}

# parameters
pars <- data.frame(chunk_index = seq_along(origin_chunks))

# submit with submit=FALSE to patch first
sjob_ttm <- slurm_apply(
  run_ttm_chunk,
  pars,
  jobname        = "ttm_full",
  nodes          = n_chunks,
  cpus_per_node  = 1,
  submit         = FALSE,
  global_objects = c("origin_chunks", "destinations"),
  pkgs           = c("r5r", "data.table"),
  slurm_options  = list(
    time          = "12:00:00",
    partition     = "normal",
    "mem-per-cpu" = "64G"
  )
)

# patch slurm_run.R -- add java options before r5r loads
log_path_full <- list.files(
  "/home/lshit9/TACTIC/build_distance_matrix",
  pattern    = "_rslurm_ttm_full",
  full.names = TRUE
)

run_script <- readLines(file.path(log_path_full, "slurm_run.R"))
writeLines(c('options(java.parameters = "-Xmx48G")', run_script),
           file.path(log_path_full, "slurm_run.R"))

# patch submit.sh -- remove --vanilla
submit_script <- readLines(file.path(log_path_full, "submit.sh"))
submit_script <- gsub("--vanilla", "", submit_script)
writeLines(submit_script, file.path(log_path_full, "submit.sh"))

# confirm patches
cat("\n--- slurm_run.R first 3 lines ---\n")
cat(head(readLines(file.path(log_path_full, "slurm_run.R")), 3), sep = "\n")
cat("\n--- submit.sh ---\n")
cat(readLines(file.path(log_path_full, "submit.sh")), sep = "\n")

# save metadata for monitor script
saveRDS(sjob_ttm,      file.path(output_path, "sjob_ttm_full.rds"))
saveRDS(log_path_full, file.path(output_path, "log_path_full.rds"))
saveRDS(n_chunks,      file.path(output_path, "n_chunks.rds"))

# submit
cat("\nsubmitting", n_chunks, "jobs...\n")
system(paste("cd", log_path_full, "&& sbatch submit.sh"))

