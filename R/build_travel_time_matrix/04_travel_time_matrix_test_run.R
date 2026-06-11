# 04_travel_time_matrix_test_run.R
# computes car travel time matrix from lsoa centroids to hospital sites
# uses r5r with england osm network
# requires england-latest.osm.pbf in Data/travel_time/routing/

library(r5r)
library(dplyr)
library(purrr)
library(data.table)

dir.create("/home/lshit9/TACTIC/build_distance_matrix/Data/routing/")
r5r_network <- build_network(data_path = "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/routing")

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



system(paste("ls -lh", routing_base))

system("df -h /tmp")

cancel_jobs()
system(paste("ls -lh", routing_base))
system("df -h /tmp")





global_objects <- c("origin_chunks", "destinations")

pars_test <- data.frame(chunk_index = 1:3)

# patch the new job folder the same way as before
sjob_test4 <- slurm_apply(
  run_ttm_chunk,
  pars_test,
  jobname        = "ttm_test4",
  nodes          = 3,
  cpus_per_node  = 1,
  submit         = FALSE,   # don't submit yet -- patch first
  global_objects = global_objects,
  pkgs           = c("r5r", "data.table"),
  slurm_options  = list(
    time          = "00:30:00",
    partition     = "normal",
    "mem-per-cpu" = "32G",
    export        = "ALL,R_JAVA_PARAMETERS=-Xmx24G"
  )
)

# find the new job folder
log_path4 <- list.files("/home/lshit9/TACTIC/build_distance_matrix",
                        pattern = "_rslurm_ttm_test4",
                        full.names = TRUE)

# patch slurm_run.R to add java option at top
run_script <- readLines



# patch slurm_run.R to add java option at top
run_script <- readLines(file.path(log_path4, "slurm_run.R"))
run_script_patched <- c('options(java.parameters = "-Xmx24G")', run_script)
writeLines(run_script_patched, file.path(log_path4, "slurm_run.R"))

# patch submit.sh to remove --vanilla
submit_script <- readLines(file.path(log_path4, "submit.sh"))
submit_script <- gsub("--vanilla", "", submit_script)
writeLines(submit_script, file.path(log_path4, "submit.sh"))

# confirm patches look right
cat("--- slurm_run.R first 5 lines ---\n")
cat(head(readLines(file.path(log_path4, "slurm_run.R")), 5), sep = "\n")
cat("\n--- submit.sh ---\n")
cat(readLines(file.path(log_path4, "submit.sh")), sep = "\n")

system(paste("cd", log_path4, "&& sbatch submit.sh"))

check_queue()
check_progress()















