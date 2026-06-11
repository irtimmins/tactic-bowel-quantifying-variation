
# collect all benchmark results
bench_files <- list.files(output_path, pattern = "^bench_", full.names = TRUE)
cat("benchmark files written:", length(bench_files), "of", nrow(pars_bench), "\n")

bench_results <- rbindlist(lapply(bench_files, fread))
bench_results <- bench_results[order(n_origins, max_duration)]

cat("\n--- benchmark results ---\n")
print(bench_results)

# summary -- mean time by n_origins and max_duration
cat("\n--- mean routing time by n_origins and max_duration ---\n")
bench_results %>%
  group_by(n_origins, max_duration) %>%
  summarise(
    mean_route_secs  = round(mean(t_route_secs), 1),
    mean_total_secs  = round(mean(t_total_secs), 1),
    mean_build_secs  = round(mean(t_build_secs), 1),
    .groups = "drop"
  ) %>%
  print()

# save summary
fwrite(bench_results, file.path(output_path, "benchmark_results.csv"))
cat("\nsaved to", file.path(output_path, "benchmark_results.csv"), "\n")

#source("/home/lshit9/TACTIC/build_distance_matrix/R/build_travel_time_matrix/run_slurm_benchmarking.R")


