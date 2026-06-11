##############################




library(readr)

ttm <- read_csv("/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/output/nearest_hospital_per_lsoa.csv")
head(ttm)
#test <- read_csv("/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/output/ttm_chunk_1.csv")


ttm_lu <- read_dta("/home/lshit9/TACTIC/build_distance_matrix/temp_travel_data/bowel_pairwise_distance_matrix_lu.dta")
head(ttm_lu)

library(dplyr)
library(ggplot2)

# get nearest hospital per lsoa from lu matrix
lu_nearest <- ttm_lu %>%
  group_by(lsoa11_code) %>%
  slice_min(total_drive_time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(nearest_site_lu = sitecode,
         travel_time_lu  = total_drive_time,
         distance_lu     = total_length)

cat("lsoas in your matrix:", nrow(ttm), "\n")
cat("lsoas in lu matrix:  ", nrow(lu_nearest), "\n")

# join on lsoa
comparison <- ttm %>%
  inner_join(lu_nearest, by = "lsoa11_code")

cat("lsoas in both:       ", nrow(comparison), "\n")

# lsoas only in yours
only_yours <- setdiff(ttm$lsoa11_code, lu_nearest$lsoa11_code)
cat("only in yours:       ", length(only_yours), "\n")

# lsoas only in lu
only_lu <- setdiff(lu_nearest$lsoa11_code, ttm$lsoa11_code)
cat("only in lu:          ", length(only_lu), "\n")

# travel time comparison
cat("\n--- travel time comparison (mins) ---\n")
cat("your matrix:\n")
print(summary(comparison$travel_time_min))
cat("lu matrix:\n")
print(summary(comparison$travel_time_lu))

# difference
comparison <- comparison %>%
  mutate(time_diff = travel_time_min - travel_time_lu)

cat("\n--- difference (yours minus lu) ---\n")
print(summary(comparison$time_diff))

# same nearest hospital?
same_site <- comparison %>%
  mutate(same = nearest_site == nearest_site_lu) %>%
  count(same) %>%
  mutate(pct = round(100 * n / sum(n), 1))

cat("\n--- same nearest hospital ---\n")
print(same_site)

# correlation
cat("\ncorrelation:", round(cor(comparison$travel_time_min,
                               comparison$travel_time_lu), 4), "\n")

# save all comparison plots to one pdf
pdf(file.path(output_path, "comparison_plots.pdf"), width = 10, height = 8)

# scatter plot
ggplot(comparison %>% sample_n(min(5000, nrow(comparison))),
       aes(x = travel_time_lu, y = travel_time_min)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  labs(x        = "lu matrix -- nearest hospital (mins)",
       y        = "your matrix -- nearest hospital (mins)",
       title    = "nearest hospital travel time comparison",
       subtitle = paste0(nrow(comparison), " lsoas")) +
  theme_minimal()

# histogram of differences
ggplot(comparison, aes(x = time_diff)) +
  geom_histogram(binwidth = 1, fill = "steelblue", colour = "white") +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
  labs(x     = "difference in travel time (yours minus lu, mins)",
       y     = "count",
       title = "distribution of travel time differences") +
  theme_minimal()

# distribution of nearest hospital travel times
ggplot(comparison %>%
         tidyr::pivot_longer(cols      = c(travel_time_min, travel_time_lu),
                             names_to  = "source",
                             values_to = "travel_time") %>%
         mutate(source = if_else(source == "travel_time_min", "yours", "lu")),
       aes(x = travel_time, colour = source)) +
  geom_density() +
  labs(x     = "nearest hospital travel time (mins)",
       y     = "density",
       title = "distribution of nearest hospital travel times") +
  theme_minimal()

dev.off()

cat("saved to", file.path(output_path, "comparison_plots.pdf"), "\n")
cat("download via WinSCP from:\n")
cat(file.path(output_path, "comparison_plots.pdf"), "\n")
# where do they disagree most?
cat("\n--- largest discrepancies ---\n")
comparison %>%
  mutate(abs_diff = abs(time_diff)) %>%
  arrange(desc(abs_diff)) %>%
  select(lsoa11_code, nearest_site, nearest_site_lu,
         travel_time_min, travel_time_lu, time_diff) %>%
  head(20) %>%
  print()

##############################################################################


# investigate the 1869 lsoas missing from your matrix
cat("lsoas in lu but not yours:\n")
missing_lsoas_detail <- lu_nearest %>%
  filter(lsoa11_code %in% only_lu) %>%
  summarise(
    n             = n(),
    min_lu_time   = min(travel_time_lu),
    median_lu_time = median(travel_time_lu),
    max_lu_time   = max(travel_time_lu)
  )
print(missing_lsoas_detail)

# are they wales lsoas?
lu_nearest %>%
  filter(lsoa11_code %in% only_lu) %>%
  mutate(country = substr(lsoa11_code, 1, 1)) %>%
  count(country)
























