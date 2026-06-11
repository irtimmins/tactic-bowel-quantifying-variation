# simple_slurm_test.R
# minimal slurm test to confirm compute nodes work correctly
# no java, no r5r -- just basic r execution

library(rslurm)

output_path <- "/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/output"
dir.create(output_path, showWarnings = FALSE)

# simplest possible function -- just return node info
simple_test <- function(i) {
  paste0("hello from job ", i,
         " on node ",    Sys.info()["nodename"],
         " at ",         format(Sys.time(), "%H:%M:%S"),
         " r version: ", R.version$version.string)
}

pars_simple <- data.frame(i = 1:3)

sjob_simple <- slurm_apply(
  simple_test,
  pars_simple,
  jobname       = "simple_test",
  nodes         = 3,
  cpus_per_node = 1,
  submit        = TRUE,
  slurm_options = list(
    time          = "00:05:00",
    partition     = "normal",
    "mem-per-cpu" = "2G"
  )
)

saveRDS(sjob_simple, file.path(output_path, "sjob_simple.rds"))

cat("jobs submitted\n")
cat("monitor with:\n")
cat('system("squeue -u lshit9")\n')

# wait a moment then check queue
Sys.sleep(10)
system("squeue -u lshit9 --format='%.10i %.15j %.8T %.10M %.5D %R'")

# poll until complete
cat("waiting for jobs to complete...\n")
for (i in 1:12) {
  Sys.sleep(10)
  n_running <- system(
    "squeue -u lshit9 -n simple_test -h | wc -l",
    intern = TRUE
  )
  cat(format(Sys.time(), "%H:%M:%S"), "jobs still in queue:", n_running, "\n")
  if (as.integer(n_running) == 0) break
}

# collect results
log_simple <- list.files(
  "/home/lshit9/TACTIC/build_distance_matrix",
  pattern    = "_rslurm_simple_test",
  full.names = TRUE
)

cat("\n--- job 0 ---\n")
cat(readLines(file.path(log_simple, "slurm_0.out")), sep = "\n")
cat("\n--- job 1 ---\n")
cat(readLines(file.path(log_simple, "slurm_1.out")), sep = "\n")
cat("\n--- job 2 ---\n")
cat(readLines(file.path(log_simple, "slurm_2.out")), sep = "\n")

# collect via rslurm
results <- get_slurm_out(sjob_simple, outtype = "raw", wait = FALSE)
cat("\nresults:\n")
print(results)


#source("/home/lshit9/TACTIC/build_distance_matrix/R/build_travel_time_matrix/simple_slurm_test.R")