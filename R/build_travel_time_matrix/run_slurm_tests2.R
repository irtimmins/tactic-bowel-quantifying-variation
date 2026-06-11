#####

# minimal_r5r_test.R
# single slurm job, 5 origins x 5 destinations
# goal: confirm r5r works on compute node

library(rslurm)

output_path <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/output"

minimal_r5r <- function(i) {

  options(java.parameters = "-Xmx24G")
  library(r5r)
  library(data.table)

  routing_base <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/routing"

  # isolated working directory
  job_dir <- file.path(tempdir(), paste0("job_", i))
  dir.create(job_dir, showWarnings = FALSE)
  old_wd <- setwd(job_dir)
  on.exit(setwd(old_wd))

  # copy network locally
  local_routing <- file.path(job_dir, "routing")
  dir.create(local_routing, showWarnings = FALSE)

  cat(format(Sys.time(), "%H:%M:%S"), "copying network\n")
  file.copy(file.path(routing_base, "england-260608.osm.pbf"),
            file.path(local_routing, "england-260608.osm.pbf"))
  file.copy(file.path(routing_base, "network.dat"),
            file.path(local_routing, "network.dat"))

  cat(format(Sys.time(), "%H:%M:%S"), "building network\n")
  r5r_network <- build_network(data_path = local_routing, verbose = FALSE)

  # tiny test -- 5 origins, 5 destinations hardcoded
  origins_test <- data.frame(
    id  = c("o1", "o2", "o3", "o4", "o5"),
    lon = c(-0.1278, -0.1000, -0.0800, -0.1500, -0.1200),
    lat = c(51.5074,  51.5100,  51.4900,  51.5200,  51.4800)
  )

  destinations_test <- data.frame(
    id  = c("d1", "d2", "d3", "d4", "d5"),
    lon = c(-0.1000, -0.0900, -0.1100, -0.1300, -0.0700),
    lat = c(51.5200,  51.5000,  51.4700,  51.5100,  51.4900)
  )

  cat(format(Sys.time(), "%H:%M:%S"), "running matrix\n")
  ttm <- travel_time_matrix(
    r5r_network,
    origins            = origins_test,
    destinations       = destinations_test,
    mode               = "CAR",
    departure_datetime = as.POSIXct("2024-06-11 10:00:00"),
    max_trip_duration  = 60L
  )

  stop_r5(r5r_network)
  rJava::.jgc(R.gc = TRUE)
  unlink(job_dir, recursive = TRUE)

  cat(format(Sys.time(), "%H:%M:%S"), "done -- rows:", nrow(ttm), "\n")
  return(ttm)
}

pars_minimal <- data.frame(i = 1)

sjob_minimal <- slurm_apply(
  minimal_r5r,
  pars_minimal,
  jobname       = "r5r_minimal",
  nodes         = 1,
  cpus_per_node = 1,
  submit        = FALSE,
  slurm_options = list(
    time          = "00:30:00",
    partition     = "normal",
    "mem-per-cpu" = "32G"
  )
)

# patch
log_minimal <- list.files(
  "/home/lshit9/TACTIC/build_distance_matrix",
  pattern    = "_rslurm_r5r_minimal",
  full.names = TRUE
)

run_script <- readLines(file.path(log_minimal, "slurm_run.R"))
writeLines(c('options(java.parameters = "-Xmx24G")', run_script),
           file.path(log_minimal, "slurm_run.R"))

submit_script <- readLines(file.path(log_minimal, "submit.sh"))
submit_script <- gsub("--vanilla", "", submit_script)
writeLines(submit_script, file.path(log_minimal, "submit.sh"))

cat("--- slurm_run.R first 3 lines ---\n")
cat(head(readLines(file.path(log_minimal, "slurm_run.R")), 3), sep = "\n")
cat("\n--- submit.sh ---\n")
cat(readLines(file.path(log_minimal, "submit.sh")), sep = "\n")

# submit
system(paste("cd", log_minimal, "&& sbatch submit.sh"))

# poll until done
cat("waiting...\n")
for (i in 1:24) {
  Sys.sleep(15)
  n_running <- as.integer(system(
    "squeue -u lshit9 -n r5r_minimal -h | wc -l",
    intern = TRUE
  ))
  cat(format(Sys.time(), "%H:%M:%S"), "jobs in queue:", n_running, "\n")
  if (n_running == 0) break
}

# read log
cat("\n--- output ---\n")
cat(readLines(file.path(log_minimal, "slurm_0.out")), sep = "\n")

# collect result
result <- get_slurm_out(sjob_minimal, outtype = "raw", wait = FALSE)
cat("\nresult:\n")
print(result)


