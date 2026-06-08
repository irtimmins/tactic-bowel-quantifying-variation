# Waiting times by CCI group - no function dependencies

wt_vars   <- c("wt_dx_to_tx", "wt_dx_to_dtt", "wt_dtt_to_tx")
wt_labels <- c("Dx -> Treatment", "Dx -> DTT", "DTT -> Treatment")

df_cci <- colon_cohort %>%
  filter(
    !is.na(wt_dx_to_tx), !is.na(wt_dx_to_dtt), !is.na(wt_dtt_to_tx),
    wt_dx_to_tx >= 0, wt_dx_to_dtt >= 0, wt_dtt_to_tx >= 0,
    !is.na(cci_group)
  )

df_cci %>%
  group_by(cci_group) %>%
  summarise(
    n               = n(),
    # Dx to Treatment
    median_dx_tx    = median(wt_dx_to_tx,  na.rm = TRUE),
    q25_dx_tx       = quantile(wt_dx_to_tx,  0.25, na.rm = TRUE),
    q75_dx_tx       = quantile(wt_dx_to_tx,  0.75, na.rm = TRUE),
    mean_dx_tx      = mean(wt_dx_to_tx,    na.rm = TRUE),
    # Dx to DTT
    median_dx_dtt   = median(wt_dx_to_dtt, na.rm = TRUE),
    q25_dx_dtt      = quantile(wt_dx_to_dtt, 0.25, na.rm = TRUE),
    q75_dx_dtt      = quantile(wt_dx_to_dtt, 0.75, na.rm = TRUE),
    mean_dx_dtt     = mean(wt_dx_to_dtt,   na.rm = TRUE),
    # DTT to Treatment
    median_dtt_tx   = median(wt_dtt_to_tx, na.rm = TRUE),
    q25_dtt_tx      = quantile(wt_dtt_to_tx, 0.25, na.rm = TRUE),
    q75_dtt_tx      = quantile(wt_dtt_to_tx, 0.75, na.rm = TRUE),
    mean_dtt_tx     = mean(wt_dtt_to_tx,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(digits = 1)
