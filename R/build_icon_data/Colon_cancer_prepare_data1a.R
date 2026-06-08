
# -----------------------------------------------------------------------------
# 1.  NCRAS: restrict to colon cancer only (C18.x), 2015+
# -----------------------------------------------------------------------------

library(arrow)
library(haven)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(purrr)

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

# Colon ICD10 (C18.x)
colon_icd <- c("C18","C180","C181","C182","C183","C184","C185","C186","C187","C188","C189")

ncras_colon <- read_parquet(
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
  filter(sitestr %in% colon_icd) %>%
  
# ncras_colon %>%
#   filter(stage_best != 4) %>%
#   group_by(ydiag) %>%
#   summarise(n = n())

# summary(as.factor(ncras_colon$stage_best))
# summary(as.factor(ncras_colon$ydiag))
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

saveRDS(
  ncras_colon,
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/ncras_colon_2015_2022.rds"
)


ncras_colon <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/ncras_colon_2015_2022.rds")

summary(as.factor(ncras_colon$route_combined))


ncras_colon_ids <- ncras_colon %>%
  distinct(pseudo_patientid) %>%
  pull(pseudo_patientid) %>%
  as.character()

cat("Colon cancer patients (NCRAS):", nrow(ncras_colon),
    "| Patients:", n_distinct(ncras_colon$pseudo_patientid), "\n")


# 2. OPCS CODE DEFINITIONS 
# ---------------------------------------------------------------------

# 4-char OPCS codes for colon major resections 
opcs_colon_codes <- list(
  right_hemi = c(
    "H061","H062","H063","H064","H065","H066","H068","H069",
    "H071","H072","H073","H074","H075","H076","H078","H079",
    "H112","H116","H118","H119"
  ),
  transverse = c("H081","H082","H083","H084","H085","H086","H088","H089"),
  left_hemi  = c("H091","H092","H093","H094","H095","H096","H098","H099","H111"),
  sigmoid    = c("H101","H102","H103","H104","H105","H106","H108","H109"),
  total_subtotal = c(
    "H041","H042","H043","H044","H045","H046","H048","H049",
    "H051","H052","H053","H054","H055","H056","H058","H059",
    "H291","H292","H293","H294","H295","H296","H298","H299",
    "H113","H114","H414"
  )
)

# Include 3-char "family" codes too (some feeds can contain them)
opcs_colon_3char <- c("H06","H07","H08","H09","H10","H04","H05","H29")

# Collapsed vectors
opcs_colon_4char <- unique(unlist(opcs_colon_codes))
opcs_colon_all   <- unique(c(opcs_colon_4char, opcs_colon_3char))

# Simple lookup: OPCS -> procedure label (only for 4-char codes)
opcs_colon_lookup <- setNames(
  rep(names(opcs_colon_codes), lengths(opcs_colon_codes)),
  unlist(opcs_colon_codes)
)

# Emergency admission method codes (HES)
admimeth_emerg <- c("21","22","23","24","25","28","2A","2B","2C","2D")

# ---------------------------------------------------------------------
# 3) Read HES APC for those patients (operations only)
# ---------------------------------------------------------------------

hes_apc_file_list <- list.files(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/HES/APC/",
  pattern    = "FILE*",
  full.names = TRUE
) %>%
  keep(~{
    yr <- str_extract(.x, "(?<=HES_APC_)\\d{4}") %>% as.integer()
    !is.na(yr) && yr %in% 2014:2024
  })

stopifnot(length(hes_apc_file_list) > 0)

op_cols      <- paste0("OPERTN_", str_pad(1:24, 2, pad = "0"))
opdate_cols  <- paste0("OPDATE_", str_pad(1:24, 2, pad = "0"))
diag_4_cols  <- paste0("DIAG_4_", str_pad(1:20, 2, pad = "0"))  # 4-char ICD, use for Charlson

hes_cols_select <- c(
  "STUDY_ID", "ADMIDATE", "ADMIMETH", "PROCODE3", "SITETRET",
  "EPISTART", "EPIORDER", "EPITYPE",
  op_cols, opdate_cols, diag_4_cols
)

hes_apc_file_list
#test <-  read_parquet("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/HES/APC/FILE0224582_NIC656757_HES_APC_202399.parquet")
#names(test)

hes_apc_raw <- map_dfr(
  hes_apc_file_list,
  ~{
    read_parquet(.x, col_select = any_of(hes_cols_select)) %>%
      filter(STUDY_ID %in% ncras_colon_ids) %>%
      mutate(
        STUDY_ID = as.character(STUDY_ID),
        ADMIMETH = as.character(ADMIMETH),
        EPISTART = as.Date(EPISTART),
        ADMIDATE = as.Date(ADMIDATE),
        across(any_of(op_cols),     as.character),
        across(any_of(opdate_cols), as.Date),
        across(any_of(diag_4_cols), as.character)
      )
  },
  .progress = TRUE
)

saveRDS(
  hes_apc_raw,
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/hes_apc_raw_colon_2014_2022.rds"
)

hes_apc_raw <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/hes_apc_raw_colon_2014_2022.rds")

# names(hes_apc_raw)
# test <- hes_apc_raw %>%
#   filter(substr(OPERTN_01, 1,1) == "H") %>%
#     slice(1:100000)
# summary(as.factor(test$OPERTN_01))
# sum(test$OPERTN_01 == "H07")
# 
# test2 <- test %>%
#   filter(substr(OPERTN_01, 1,3) == "H07")
# summary(as.factor(test2$OPERTN_01))

cat("HES rows (restricted to colon NCRAS IDs):", nrow(hes_apc_raw), "\n")

#op_cols <- names(hes_apc_raw)[startsWith(names(hes_apc_raw), "OPERTN_")]
#op_cols

# ---------------------------------------------------------------------
# 4. Pivot operations long and keep only colon cancer resection OPCS
# ---------------------------------------------------------------------
# sum(substr(hes_apc_raw$SITETRET, 1,1) == "R")
# sum(substr(hes_apc_raw$PROCODE3,1,1) == "R")
# sum(grepl("^R[A-Z0-9]{4}$", hes_apc_raw$SITETRET))
# 
# filter1 <- substr(hes_apc_raw$SITETRET, 1,1) == "R"
# filter2 <- grepl("^R[A-Z0-9]{4}$", hes_apc_raw$SITETRET)
# test3 <- hes_apc_raw[filter1 & !(filter2),]
# summary(as.factor(test3$SITETRET))
# summary(as.factor(test3$PROCODE3))
#summary(as.factor(substr(hes_apc_raw$PROCODE3, 1,3)))
#length(unique(hes_apc_raw$PROCODE3[substr(hes_apc_raw$PROCODE3, 1,1) == "R"]))


names(hes_apc_raw)

hes_opcs_long <- 
  hes_apc_raw %>%
  select(any_of(hes_cols_select)) %>%
  filter(grepl("^R[A-Z0-9]{4}$", SITETRET)) %>%
  filter(substr(PROCODE3, 1,1) == "R") %>%
  filter(EPITYPE == "1", !is.na(OPERTN_01), OPERTN_01 != "-") %>%
 # slice(1:10000) %>% View()
  pivot_longer(
    cols      = all_of(op_cols),
    names_to  = "op_position",
    values_to = "opcs_code"
  ) %>%
  filter(!is.na(opcs_code), opcs_code != "-") %>%
  #   slice(1:10000) %>%
  # View()
  mutate(
    opcs3 = str_sub(opcs_code, 1, 3),
    op_position_n = as.integer(str_extract(op_position, "[0-9]+")),
    colon_proc_type = case_when(
      opcs_code %in% opcs_colon_4char ~ unname(opcs_colon_lookup[opcs_code]),
      opcs3 %in% opcs_colon_3char     ~ opcs3,  # keep family code if it appears
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(colon_proc_type)) %>%
  mutate(
    # if 3-char slipped through, map it to groups for consistency
    colon_proc_type = case_when(
      colon_proc_type %in% c("H06","H07") ~ "right_hemi",
      colon_proc_type == "H08"            ~ "transverse",
      colon_proc_type == "H09"            ~ "left_hemi",
      colon_proc_type == "H10"            ~ "sigmoid",
      colon_proc_type %in% c("H04","H05","H29") ~ "total_subtotal",
      TRUE ~ colon_proc_type
    )
  )

#summary(hes_opcs_long$EPISTART)

# ---------------------------------------------------------------------
# 5. Collapse to episode-level colon major resections
# ---------------------------------------------------------------------

first_or_na <- function(x) if (length(x) == 0) NA_character_ else dplyr::first(x)

# Optional: stable primary type priority (if multiple colon codes in an episode)
proc_priority <- c("total_subtotal", "right_hemi", "left_hemi", "sigmoid", "transverse")

#dplyr::first(c("a", "b", "c"))

names(hes_opcs_long)

# hes_opcs_long %>%
#   group_by(op_position) %>%
#   summarise(n = n())


hes_colon_episodes <-
  hes_opcs_long %>%
  select(!starts_with("DIAG")) %>%
  arrange(STUDY_ID, EPISTART, op_position_n) %>%
  group_by(STUDY_ID, ADMIDATE, EPISTART, EPIORDER, EPITYPE, PROCODE3, SITETRET, ADMIMETH)  %>%
  mutate(primary_flag = op_position_n == min(op_position_n)) %>%
 # filter(n() > 1) %>% View()
# summary(as.factor(hes_colon_episodes$op_position_n))
# length(unique(hes_colon_episodes$STUDY_ID))
# hes_colon_episodes %>% View()
  # summarise(
  #   # Choose a primary type with a clear hierarchy
  #   colon_proc_type = first_or_na(intersect(proc_priority, unique(colon_proc_type))),
  #   
  #   # Choose the first matched OPCS code in operation order (good proxy for "primary")
  #   colon_opcs_primary = first_or_na(opcs_code),
  #   
  #   # Keep all matched codes for QA
  #   all_colon_opcs = paste(unique(opcs_code), collapse = "; "),
  #   
  #   .groups = "drop"
  # ) %>%
  # 
  mutate(
    all_colon_opcs          = paste(unique(opcs_code),       collapse = ";"),
    all_colon_proc_types    = paste(unique(colon_proc_type), collapse = ";"),
    all_colon_opcs_fields = paste(unique(op_position),       collapse = ";"),
    n_colon_codes_in_episode = n()
  ) %>%
  # Now drop to one row: the primary (lowest op_position_n) 
  filter(primary_flag) %>%
  slice(1) %>%   # guard: if two codes share the same min position, take the first
  ungroup() %>%
  rename(
    colon_opcs_primary          = opcs_code,
    colon_proc_type_primary     = colon_proc_type,
    colon_opcs_primary_position = op_position_n
  ) %>%
  mutate(
    emergency = ADMIMETH %in% admimeth_emerg
  )
    # filter(n()>1) %>%
    # slice(1:10000) %>% View()
# hes_colon_episodes %>% filter(n_colon_codes_in_episode > 1) %>%
#   slice(1:10000) %>% View()
cat("Colon major resection episodes:", nrow(hes_colon_episodes),
    "| Patients:", n_distinct(hes_colon_episodes$STUDY_ID), "\n")

# hes_opcs_long %>%
#   filter(STUDY_ID == "1321793") %>%
#   View()
