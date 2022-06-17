# MUSIC GET CLINICAL DATA
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
output_dir <- "data/music/clinical/"
# -----------------------------------------------------------------------------
# Description
# -----------------------------------------------------------------------------
# This file queries the API to retrieve MUSIC data. Data exported into .csv and
# .rds formats.
#
# Credentials are read from 'credentials.json'. MUSIC_API_KEY is the key that
# this script looks for.
#
# Initial data read into df.raw. df_demo and df_main are created.
# Couple of transformations to squash df_main.
# Finally cleaned file loaded into df.
# -----------------------------------------------------------------------------
("Loading required packages...")
library(dplyr)
library(jsonlite)
library(httr)

# Retrieve data from API
# -----------------------------------------------------------------------------
# API Key is retrieved from file credentials.json

print("Retrieving data from RedCap API...")
music_api_key <- jsonlite::read_json(path = "credentials.json")$MUSIC_API_KEY
url <- "https://ecrf.igmm.ed.ac.uk/api/"
form_data <- list(
  "token" = music_api_key,
  content = "record",
  action = "export",
  format = "json",
  type = "flat",
  csvDelimiter = "",
  rawOrLabel = "raw",
  rawOrLabelHeaders = "raw",
  exportCheckboxLabel = "false",
  exportSurveyFields = "false",
  exportDataAccessGroups = "false",
  returnFormat = "json"
)
response <- httr::POST(url, body = form_data, encode = "form")
result <- content(response, as = "text")
df_raw <- jsonlite::fromJSON(result)

# Reformat data to cleaned dataframe
# -----------------------------------------------------------------------------
# Create 2 data frames - df_demo (with the grouping info)
# and df_main (with the repeated measurements)
# Select desired columns into df_main
# Restructure df so that repeated measures will have grouping info on each line.

df_demo <- subset(df_raw, redcap_repeat_instrument == "")
df_demo <- subset(df_demo, redcap_event_name == "timepoint_1_arm_1")
df_demo <- df_demo %>% select(
  c(
    "study_id",
    "patientdetailsage",
    "sex",
    "study_group_name",
    "date_of_diagnosis",
    "patient_diagnose_age"
  )
)

# colnames(df) will give you all data
# Select the longitudinal variables you want.
df_main <- df_raw %>%
  select(c(
    "study_id",
    "redcap_repeat_instrument",
    "redcap_event_name",
    "visit_review_date",
    "total_cucq32_scored",
    "bmi_weight",
    "clinician_perform_stat",
    "hbi_total_calculation",
    "sccai_total_calculation",
    "mayo_total_score",
    "new_flare_up",
    "unplanned_hospital",
    "haemoglobin",
    "red_cell_count",
    "white_cell_count",
    "neutrophils_lab",
    "lymphocytes_lab",
    "monocytes_lab",
    "eosinophils_lab",
    "basophils_lab",
    "plt_lab",
    "crp",
    "albumin",
    "calprotectin",
    "mayo_endoscopic_findings",
    "uceis_total_score",
    "uceis_noncalc",
    "cdeis_total_score",
    "sescd_noncalc",
    "montrlcdlocation_montrl",
    "montrlcdbehav_montrl",
    "montreal_perianal",
    "montrlucextent_montrl",
    "montrl_uc_server"
  ))

# Now compress df_main into usable format
df_main <- subset(
  df_main,
  redcap_repeat_instrument != "research_samples_obtained"
)
# Fill all blanks with NA
df_main <- df_main %>% mutate_all(na_if, "")

# Drop columns where all values are NA
df_main_split_1 <- subset(df_main, redcap_repeat_instrument != "")
df_main_split_1 <- df_main_split_1[colSums(!is.na(df_main_split_1)) > 0]
df_main_split_2 <- subset(df_main, is.na(redcap_repeat_instrument))
df_main_split_2 <- df_main_split_2[colSums(!is.na(df_main_split_2)) > 0]

# Merge and squash data frame
df_main_squashed <- full_join(
  df_main_split_1,
  df_main_split_2,
  by = c("study_id", "redcap_event_name")
)

# Finally, return the cleaned data frame as df

df <- left_join(df_main_squashed, df_demo, by = "study_id")
df <- df %>% select(-c("redcap_repeat_instrument"))

# Convert factors
# -----------------------------------------------------------------------------
print("Factor conversions...")
df$redcap_event_name <- as.factor(df$redcap_event_name)
df$redcap_event_name <- recode(
  df$redcap_event_name,
  "timepoint_1_arm_1" = "timepoint_1",
  "timepoint_2_arm_1" = "timepoint_2",
  "timepoint_3_arm_1" = "timepoint_3",
  "timepoint_4_arm_1" = "timepoint_4",
  "timepoint_5_arm_1" = "timepoint_5"
)

df$sex <- as.factor(df$sex)
df$sex <- recode(df$sex, "1" = "male", "2" = "female")

df$study_group_name <- as.factor(df$study_group_name)
df$study_group_name <- recode(df$study_group_name, "1" = "CD", "2" = "UC")

df$clinician_perform_stat <- as.factor(df$clinician_perform_stat)
df$clinician_perform_stat <- recode(
  df$clinician_perform_stat,
  "1" = "Remission", "2" = "Mild", "3" = "Moderate", "4" = "Severe"
)

df$new_flare_up <- as.factor(df$new_flare_up)
df$new_flare_up <- recode(
  df$new_flare_up,
  "1" = "Yes", "0" = "No", "3" = "Ongoing"
)

df$unplanned_hospital <- as.factor(df$unplanned_hospital)
df$unplanned_hospital <- recode(df$unplanned_hospital, "1" = "Yes", "0" = "No")

df$montrlcdlocation_montrl <- as.factor(df$montrlcdlocation_montrl)
df$montrlcdlocation_montrl <- recode(df$montrlcdlocation_montrl,
                                     "0" = "L1",
                                     "1" = "L2",
                                     "2" = "L3",
                                     "3" = "L4",
                                     "4" = "L1+L4",
                                     "5" = "L2+L4",
                                     "6" = "L3+L4"
                                     )

df$montrlcdbehav_montrl <- as.factor(df$montrlcdbehav_montrl)
df$montrlcdbehav_montrl <- recode(df$montrlcdbehav_montrl, "0" = "B1", "1" = "B2", "2" = "B3", "3" = "p")

df$montreal_perianal <- as.factor(df$montreal_perianal)
df$montreal_perianal <- recode(df$montreal_perianal, "1" = "Perianal disease")

df$montrlucextent_montrl <- as.factor(df$montrlucextent_montrl)
df$montrlucextent_montrl <- recode(df$montrlucextent_montrl, "0" = "E1", "1" = "E2", "2" = "E3")

df$montrl_uc_server <- as.factor(df$montrl_uc_server)
df$montrl_uc_server <- recode(df$montrl_uc_server,"0" = "S0", "1" = "S1", "2" = "S2", "3" = "S3")

# Convert numerics
# -----------------------------------------------------------------------------
print("Numeric conversions...")
numeric_cols <- c(
  "haemoglobin",
  "red_cell_count",
  "white_cell_count",
  "neutrophils_lab",
  "lymphocytes_lab",
  "monocytes_lab",
  "basophils_lab",
  "eosinophils_lab",
  "plt_lab",
  "crp",
  "albumin",
  "total_cucq32_scored",
  "bmi_weight",
  "hbi_total_calculation",
  "sccai_total_calculation",
  "mayo_total_score",
  "uceis_total_score",
  "uceis_noncalc",
  "cdeis_total_score",
  "sescd_noncalc",
  "patientdetailsage",
  "patient_diagnose_age",
  "mayo_endoscopic_findings"
)
df[, numeric_cols] <- lapply(df[, numeric_cols], as.numeric)

# Calprotectin Cleaning
clean_calprotectin <- function(x) {
  if (is.na(x)) {
    cleaned_value <- NA
  } else {
    if (substr(x, 1, 1) == "<" | substr(x, 1, 1) == ">") {
      cleaned_value <- substring(x, 2, nchar(x))
    } else {
      cleaned_value <- x
    }
  }

  return(cleaned_value)
}
df$calprotectin <- lapply(df$calprotectin, clean_calprotectin)
df$calprotectin <- as.numeric(df$calprotectin)

# Convert dates
# -----------------------------------------------------------------------------
print("Date conversions...")
date_cols <- c("date_of_diagnosis", "visit_review_date")
df[, date_cols] <- lapply(df[, date_cols], as.Date)

# Create location column
# -----------------------------------------------------------------------------
add_center_name <- function(x) {
  if (is.na(x)) {
    cleaned_value <- NA
  } else {
    if (substr(x, 1, 2) == "90") {
      cleaned_value <- "Edinburgh"
    } else if (substr(x, 1, 2) == "91") {
      cleaned_value <- "Glasgow"
    } else if (substr(x, 1, 2) == "92") {
      cleaned_value <- "Dundee"
    } else {
      cleaned_value <- "Error"
    }
  }
  return(cleaned_value)
}

df$center <- lapply(df$study_id, add_center_name)
df$center <- as.character(df$center)
# Rearrange data frame
# -----------------------------------------------------------------------------

df <- df %>% relocate(center, .before = study_id)
df <- df %>% relocate(study_group_name, .after = study_id)
df <- df %>% relocate(visit_review_date, .after = redcap_event_name)
df <- df %>% relocate(patientdetailsage, .after = study_group_name)
df <- df %>% relocate(sex, .after = patientdetailsage)

# Rename columns
# -----------------------------------------------------------------------------

df <- df %>% dplyr::rename(age = patientdetailsage,
  cucq_total = total_cucq32_scored,
  sccai_total = sccai_total_calculation,
  hbi_total = hbi_total_calculation,
  mayo_total = mayo_total_score,
  physician_global_assessment = clinician_perform_stat,
  montreal_cd_location = montrlcdlocation_montrl,
  montreal_cd_behaviour = montrlcdbehav_montrl,
  montreal_uc_extent = montrlucextent_montrl,
  montreal_uc_severity = montrl_uc_server
)

# Save RDS and CSV files.
# Use RDS if for further R analysis to preserve data types.
# -----------------------------------------------------------------------------
current_date <- Sys.Date()
saveRDS(df, file = paste0(output_dir, "music_clinical_", current_date, ".rds"))
write.csv(
  df,
  file = paste0(output_dir, "music_clinical_", current_date, ".csv"),
  row.names = FALSE
)
print("Saving complete.")
print("Script completed successfully.")
