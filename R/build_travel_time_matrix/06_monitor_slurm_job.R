# 06 
# monitors slurm job progress and collects results when complete
# can be sourced repeatedly to check progress
# requires 04_slurm_matrix_full.R to have been run first

library(data.table)
library(dplyr)
library(haven)

base_path   <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data"
output_path <- file.path(base_path, "output")
rds_path    <- file.path(base_path, "rds")

# load run metadata
n_chunks     <- readRDS(file.path(output_path, "n_chunks.rds"))
log_path_full <- readRDS(file.path(output_path, "log_path_full.rds"))

if (!exists("origins")) {
  origins <- readRDS(file.path(rds_path, "origins.rds"))
}

# queue status
check_queue <- function() {
  system("squeue -u lshit9 -n ttm_full --format='%.10i %.15j %.8T %.10M %.5D %R'")
}

# sacct summary
check_jobs <- function() {
  system(paste0("sacct -S ", Sys.Date(),
                " -u lshit9 -j ttm_full",
                " --format=JobID,Jobname,partition,state,elapsed,ncpus -X"))
}

# chunk file progress
# fixed elapsed time check -- -n flag not supported in all sacct versions
check_progress <- function() {
  chunk_files <- list.files(output_path, pattern = "^ttm_chunk_", full.names = FALSE)
  completed   <- as.integer(gsub("ttm_chunk_|\\.csv", "", chunk_files))
  all_chunks  <- seq_len(n_chunks)
  missing     <- setdiff(all_chunks, completed)

  cat("chunk files written:", length(chunk_files), "of", n_chunks, "\n")
  cat("pct complete:       ", round(100 * length(chunk_files) / n_chunks, 1), "%\n")

  if (length(missing) > 0) {
    cat("missing chunks:     ", length(missing), "\n")
    cat("first few missing:  ", head(missing), "\n")
  } else {
    cat("all chunks complete\n")
  }
}

# failed jobs
check_failed <- function() {
  failed <- system(
    paste0("sacct -S ", Sys.Date(),
           " -u lshit9 -n ttm_full",
           " --format=JobID,state,elapsed --noheader -X",
           " | grep -i failed"),
    intern = TRUE
  )
  if (length(failed) == 0) {
    cat("no failed jobs\n")
  } else {
    cat("failed jobs:\n")
    cat(failed, sep = "\n")
  }
}

# read a specific job log
read_log <- function(i) {
  cat(readLines(file.path(log_path_full, paste0("slurm_", i, ".out"))), sep = "\n")
}

# collect results once all chunks complete
collect_results <- function() {
  chunk_files <- list.files(output_path, pattern = "^ttm_chunk_",
                            full.names = TRUE)

  if (length(chunk_files) < n_chunks) {
    cat("only", length(chunk_files), "of", n_chunks, "chunks available\n")
    cat("collect anyway? run collect_results(force=TRUE)\n")
    return(invisible(NULL))
  }

  cat("reading", length(chunk_files), "chunk files\n")
  ttm <- rbindlist(lapply(chunk_files, fread))

  cat("total matrix rows:    ", nrow(ttm), "\n")
  cat("unique origins:       ", n_distinct(ttm$from_id), "\n")
  cat("unique destinations:  ", n_distinct(ttm$to_id), "\n")
  cat("travel time range:    ", range(ttm$travel_time_p50), "\n")

  # lsoas with no reachable hospital
  missing_lsoas <- setdiff(origins$id, unique(ttm$from_id))
  cat("lsoas with no result: ", length(missing_lsoas), "\n")

  # nearest hospital per lsoa
  nearest <- ttm %>%
    group_by(from_id) %>%
    slice_min(travel_time_p50, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(lsoa11_code     = from_id,
           nearest_site    = to_id,
           travel_time_min = travel_time_p50)

  cat("lsoas with nearest:   ", nrow(nearest), "\n")

  # save
  fwrite(ttm,     file.path(output_path, "ttm_lsoa_hospital.csv"))
  fwrite(nearest, file.path(output_path, "nearest_hospital_per_lsoa.csv"))

  cat("saved:\n")
  cat(" ", file.path(output_path, "ttm_lsoa_hospital.csv"), "\n")
  cat(" ", file.path(output_path, "nearest_hospital_per_lsoa.csv"), "\n")

  return(invisible(list(ttm = ttm, nearest = nearest)))
}

# cancel all jobs if needed
cancel_jobs <- function() {
  system("scancel -u lshit9")
  Sys.sleep(5)
  n <- as.integer(system("squeue -u lshit9 -h | wc -l", intern = TRUE))
  cat("jobs remaining in queue:", n, "\n")
}

# print current status on source
cat("=== travel time matrix job monitor ===\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
check_queue()
cat("\n")
check_progress()
cat("\nfunctions available:\n")
cat("  check_queue()    -- live queue\n")
cat("  check_jobs()     -- sacct summary\n")
cat("  check_progress() -- chunk files written\n")
cat("  check_failed()   -- any failed jobs\n")
cat("  read_log(i)      -- read log for job i\n")
cat("  collect_results()-- combine chunks into final matrix\n")
cat("  cancel_jobs()    -- cancel everything\n")