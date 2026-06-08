


ncras_cols <- c(
  # Identity and linkage
  "pseudo_patientid", "pseudo_tumourid",
  
  # Diagnosis date and year
  "diagmdy", "ydiag",
  
  # Tumour characteristics
  "cancer", "sitestr", "typestr",
  "basisofdiagnosis",
  "grade",
  
  # Staging
  "stage_best", "stage_best_system",
  "t_best", "n_best", "m_best",
  "t_path", "n_path", "m_path",
  
  # Patient characteristics
  "sex", "agediag", "birthmdy",
  "ethnicity_group_broad",
  
  # Geography and deprivation
  "lsoa11_code",
  "NHSE_reversed_imd_quintile_lsoas",
  "canalliance_2024_code", "canalliance_2024_name",
  
  # Organisation
  "diag_trust", "diag_trust_name",
  "first_trust", "first_trust_name", "first_hosp_date",
  "diag_hosp" , "diag_hosp_name",
  
  
  # Pathway and waiting times
  "route_bjc", "final_route", "route_code",
  "tww_to_treat",
  
  # Treatment flags
  "sg_flag", "rt_flag", "ct_flag",
  
  # Screening
  "screendetected",
  
  # Survival
  "dead", "finmdy", "dco",
  
  # Breast specific (kept for compatibility)
  "er_status", "pr_status", "her2_status", "laterality",
  
  # Bowel specific
  "dukes", "nodesexcised", "nodesinvolved"
)

ncras <- read_parquet(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/NCRAS/NCRAS_clean_1995_2022_route.parquet"
) %>%
  mutate(across(where(is.labelled), as_factor)) %>%
  select(all_of(ncras_cols)) %>%
  mutate(
    pseudo_patientid = as.character(pseudo_patientid),
    diagmdy          = as.Date(diagmdy),
    
    final_route_chr  = na_if(as.character(final_route), ""),
    route_bjc_chr    = na_if(as.character(route_bjc), ""),
    # precedence: final_route > route_bjc
    route_combined   = coalesce(final_route_chr, route_bjc_chr),
    route_combined   = if_else(is.na(route_combined), "Unknown", route_combined),
    route_combined   = factor(route_combined)
  ) %>%
  zap_labels() %>%
  zap_formats() %>%
  filter(ydiag >= 2015) %>%
  # Colon only
  # Stage/route exclusions
  filter(!stage_best %in% c("X", "U")) %>%
  mutate(stage = case_when(
    str_starts(stage_best, "1") ~ "1",
    str_starts(stage_best, "2") ~ "2",
    str_starts(stage_best, "3") ~ "3",
    str_starts(stage_best, "4") ~ "4",
    TRUE ~ NA_character_
  )) %>%
  filter(stage %in% c("1","2","3")) %>%
  filter(!route_combined %in% c("Emergency presentation", "Unknown"))

###############################################

ncras %>%
  group_by(cancer) %>%
  summarise(n= n()) %>%
  arrange(-n) %>%
  View()

ncras %>%
 # filter(cancer == "colon") %>%
  filter(sitestr %in% og_icd10) %>%
  group_by( route_combined, ydiag) %>%
  summarise(n = n()) %>% View()


ncras %>%
  filter(cancer == "oesophagus") %>%
  group_by( route_combined, ydiag) %>%
  summarise(n = n()) %>% View()

ncras %>%
  filter(cancer == "oesophagus") %>%
  group_by( route_bjc, ydiag) %>%
  summarise(n = n()) %>% View()


ncras %>%
  filter(sitestr %in% og_icd10) %>%
  group_by( final_route, ydiag) %>%
  summarise(n = n()) %>% View()

ncras %>%
  filter(cancer == "colon") %>%
  group_by( final_route, ydiag) %>%
  summarise(n = n()) %>% View()

ncras %>%
  filter(cancer == "stomach") %>%
  group_by(ydiag, final_route) %>%
  summarise(n = n()) %>%
  mutate(percentage = (n / sum(n)) * 100) %>%
  print(n = 50)


ncras %>%
  group_by(cancer, ydiag) %>%
  summarise(n = n()) %>%
  mutate(percentage = (n / sum(n)) * 100) %>%
  group_by(cancer, ydiag) %>%
  mutate(patients = sum(n)) %>%
  arrange(-patients) %>%
  #  filter(final_route == "TWW") %>%
  filter(cancer %in% c("stomach", "oesophagus", "rectum", "ovarian")) %>%
  View()

ncras %>%
  group_by(cancer, ydiag, final_route) %>%
  summarise(n = n()) %>%
  mutate(percentage = (n / sum(n)) * 100) %>%
  filter(cancer %in% c("colon")) %>%
  View()
  group_by(cancer, ydiag) %>%
  mutate(patients = sum(n)) %>%
  arrange(-patients) %>%
  filter(final_route == "TWW") %>%
  filter(cancer %in% c("colon")) %>%
  #View()
 # filter(patients > 1000) %>%
  ggplot(aes(x = ydiag, y = percentage, colour = cancer))+
  theme_classic()+
  geom_line()+
  scale_x_continuous(breaks = 2015:2022)+
  scale_y_continuous(limits = c(0, NA))


################################################

colon_cohort <- readRDS(
"E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_2015_2022.rds"
)

df <- colon_cohort %>%
  filter(
    !is.na(wt_dx_to_tx), !is.na(wt_dx_to_dtt), !is.na(wt_dtt_to_tx),
    wt_dx_to_tx  > 0, wt_dx_to_tx <= 180,
    wt_dx_to_dtt > 0,
    wt_dtt_to_tx > 0
  ) %>%
  mutate(
    ydiag     = as.integer(ydiag),
    age_group = cut(agediag, breaks = c(0, 49, 59, 69, 79, Inf),
                    labels = c("<50", "50-59", "60-69", "70-79", "80+"),
                    right = TRUE),
    season = case_when(
      month(diagmdy) %in% c(12, 1, 2)  ~ "Winter",
      month(diagmdy) %in% c(3, 4, 5)   ~ "Spring",
      month(diagmdy) %in% c(6, 7, 8)   ~ "Summer",
      month(diagmdy) %in% c(9, 10, 11) ~ "Autumn"
    ),
    season    = factor(season, levels = c("Spring", "Summer", "Autumn", "Winter")),
    diag_hosp = factor(diag_hosp)
  ) %>%
  filter(!is.na(wt_dx_to_dtt), !is.na(diag_hosp)) %>%
  mutate(ydiag = factor(ydiag)) %>%
  mutate(org_trust_treat_start = substr(org_treat_start, 1, 3)) %>%
  mutate(diag_surgery_hes = diag_trust == substr(PROCODE3, 1, 3),
         diag_surgery_cwt = diag_trust == org_trust_treat_start,
         diag_surgery_hes_hosp = as.character(diag_hosp) == as.character(SITETRET)) %>%
  mutate(surgery_hes_cwt = substr(PROCODE3, 1, 3) == org_trust_treat_start) %>%
  mutate(dtt_trust = substr(org_dec_to_treat, 1, 3),
         dtt_treatment_change = dtt_trust == substr(PROCODE3, 1, 3)) 

#df %>%
#  select(dtt_trust, diag_surg)
#names(df)

df %>%
  group_by(dtt_treatment_change) %>%
  summarise(mean_wt = mean(wt_dx_to_tx),
              n = n())

df %>%
  group_by(route_combined) %>%
  summarise(mean_wt = mean(wt_dx_to_tx),
            n = n())


df %>%
  group_by(diag_surgery_hes) %>%
  summarise(n = n())

df %>%
  group_by(diag_surgery_cwt) %>%
  summarise(n = n())

df %>%
  group_by(surgery_hes_cwt) %>%
  summarise(n = n())

df %>%
 # filter(nchar(as.character(diag_hosp)) == 5,
#         nchar(as.character(SITETRET)) == 5) %>%
  group_by(diag_surgery_hes_hosp) %>%
  summarise(n = n())

df$ydiag

df %>%
  mutate(ydiag = as.integer(as.character(ydiag))) %>%
  group_by(ydiag, dtt_treatment_change) %>%
    summarise(wt = mean(wt_dx_to_dtt),
            wt_sd = sd(wt_dx_to_dtt),
            n =  n()) %>%
  ungroup() %>%
 # mutate(ydiag = as.numeric(ydiag)) #%>%
  ggplot(aes(x = ydiag, y = wt, colour = diag_surgery_hes))+
  theme_classic()+
  geom_line()+
  scale_y_continuous( limits = c(0, NA))

library(ggplot2)
df %>%
  mutate(ydiag = as.integer(as.character(ydiag))) %>%
  group_by(ydiag, route_combined) %>%
  summarise(wt = mean(wt_dx_to_dtt),
            wt_sd = sd(wt_dx_to_dtt),
            n =  n()) %>%
  ungroup() %>%
  # mutate(ydiag = as.numeric(ydiag)) #%>%
  ggplot(aes(x = ydiag, y = wt, colour = route_combined))+
  theme_classic()+
  geom_line()+
  scale_y_continuous( limits = c(0, NA))


length(unique(df$SITETRET))
length(unique(df$PROCODE3))
length(unique(df$diag_hosp))

df$diag_hosp
df %>%
  group_by(ydiag, route_combined) %>%
  summarise(n = n()) %>% View()



df %>%
  group_by(diag_surgery_hes) %>%
  summarise(wt = mean(wt_dx_to_dtt),
            wt_sd = sd(wt_dx_to_dtt),
            n =  n())



summary(df$wt_dx_to_tx)
names(df)
df %>%
  group_by(ydiag, route_combined) %>%
  summarise(wt = mean(wt_dx_to_tx),
            wt_sd = sd(wt_dx_to_tx),
            n =  n()) %>%
  View()
summary(df$wt_dx_to_tx)

summary(as.factor(df$diag_trust))
summary(as.factor(df$PROCODE3))

df %>%
  mutate(change_trust = 1*(diag_trust == PROCODE3)) %>%
  group_by(ydiag, change_trust) %>%
  summarise(wt = mean(wt_dx_to_tx),
            wt_sd = sd(wt_dx_to_tx),
            n =  n())


df %>%
  mutate(change_trust = 1*(diag_trust == PROCODE3)) %>%
 # filter(ydiag %in% 2020:2022) %>%
  group_by(route_combined, change_trust) %>%
  summarise(wt = mean(wt_dx_to_tx),
            wt_sd = sd(wt_dx_to_tx),
            n =  n()) %>%
  View()



df %>%
  group_by(ydiag, NHSE_reversed_imd_quintile_lsoas) %>%
  summarise(wt = mean(wt_dx_to_tx),
            wt_sd = sd(wt_dx_to_tx),
            n =  n()) %>%
  View()

df %>%
  filter(ydiag %in% 2020:2020) %>%
  group_by(route_combined) %>%
  summarise(wt = mean(wt_dx_to_tx),
            wt_sd = sd(wt_dx_to_tx),
            n =  n())
summary(df$wt_dx_to_tx)






