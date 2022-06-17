# GIDAMPS GET CLINICAL DATA
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
output_dir <- "data/gidamps/clinical/"
# -----------------------------------------------------------------------------
# Description
# -----------------------------------------------------------------------------
# This file queries the API to retrieve MUSIC data. Data exported into .csv and 
# .rds formats.
#
# Credentials are read from 'credentials.json'. MUSIC_API_KEY is the key that
# this script looks for.
#
# Initial data read into df.raw. df.demo and df.main are created.
# Couple of transformations to squash df.main.
# Finally cleaned file loaded into df.
# -----------------------------------------------------------------------------
print("Loading required packages...")
library(dplyr)
library(jsonlite)
library(httr)
print("Complete.")

print("Retrieving data from RedCap API...")

# API Key is retrieved from file credentials.json
gidamps_api_key <- jsonlite::read_json(path="credentials.json")$GIDAMPS_API_KEY
url <- "https://ecrf.igmm.ed.ac.uk/api/"
formData <- list("token"=gidamps_api_key,
                 content='record',
                 action='export',
                 format='json',
                 type='flat',
                 csvDelimiter='',
                 'forms[0]'='gidamps_participant_questionnaire',
                 'forms[1]'='ibd_background_clinician',
                 'forms[2]'='patient_details',
                 'forms[3]'='sampling',
                 rawOrLabel='raw',
                 rawOrLabelHeaders='raw',
                 exportCheckboxLabel='false',
                 exportSurveyFields='false',
                 exportDataAccessGroups='false',
                 returnFormat='json'
)
response <- httr::POST(url, body = formData, encode = "form")
result <- content(response, as="text")
df <- jsonlite::fromJSON(result)
print("Complete.")

# Create a clean demographic data frame

# Stage 1 clean (separate sampling columns out)
print("Creating and cleaning demographic dataframe...")
demographic_df <- df %>%
  filter(redcap_repeat_instrument == "") %>%
  select(-starts_with("blood_sample_")) %>%
  select(-starts_with("faecal_sample_")) %>%
  select(-starts_with("biopsy_sample_")) %>%
  select(-starts_with("sccai_")) %>%
  select(-starts_with("hbi")) %>%
  select(-ends_with("_lab")) %>%
  select(-ends_with("_test")) %>%
  select(-ends_with("_total")) %>%
  select(-ends_with("_experiment")) %>%
  select(-starts_with("drug_level")) %>%
  select(-starts_with("sampling_")) %>%
  select(-starts_with("endoscopy")) %>%
  select(-starts_with("radiology")) %>%
  select(-starts_with("steroids")) %>%
  select(-ends_with("steroids")) %>%
  select(-starts_with("rescue_")) %>%
  select(-starts_with("mri")) %>%
  select(-starts_with("antibiotics")) %>%
  select(
    -c(
      redcap_repeat_instrument,
      redcap_repeat_instance,
      email,
      consent,
      study_group,
      study_group_hc,
      bloodtestdate_as_lbt,
      faecal_test_date,
      sampls_date_experiment4,
      patient_active_symptomyn1,
      baseline_gi_symptoms_desc,
      blood_test_date_red,
      haemoglobin,
      red_cell_count,
      white_cell_count,
      crp,
      albumin,
      faecal_test_date_2,
      ibd_status,
      radilogy_r_ctabdo_pelvic,
      ct_abdomen_and_or_pelvis,
      calprotectin,
      iv_steroids_day,
      date_test_adali_inflixi,
      medication_comments,
      histopathology_report_text,
      )
    )

as.matrix(colnames(demographic_df))

# Stage 2 clean (Rename columns.)
print("Renaming columns for demographic dataframe...")
cleaned_demographic_df <- demographic_df %>%
  rename(
    is_healthy_control = pd_gi_participant_category,
    group = study_group_name,
    consent_date = date_enrolledconsent,
    height = bmi_height,
    weight = bmi_weight,
    age_at_diagnosis = diagnosis_age,
    is_new_diagnosis = bl_new_diagnosis,
    montreal_l1 = baseline_mont_cd_loc___0,
    montreal_l2 = baseline_mont_cd_loc___1,
    montreal_l3 = baseline_mont_cd_loc___2,
    montreal_l4 = baseline_mont_cd_loc___3,
    montreal_b1 = baseline_mont_cd_beh___0,
    montreal_b2 = baseline_mont_cd_beh___1,
    montreal_b3 = baseline_mont_cd_beh___2,
    montreal_p = baseline_mont_cd_beh___3,
    montreal_extent = baseline_mont_uc_extent,
    montreal_severity = baseline_mont_seve_uc,
    eims_arthralgia = baseline_eims___2,
    eims_ank_spon = baseline_eims___3,
    eims_erythema_nodosum = baseline_eims___5,
    eims_pyoderma = baseline_eims___6,
    eims_uveitis = baseline_eims___10,
    eims_scleritis = baseline_eims___12,
    eims_sacroileitis = baseline_eims___8,
    eims_none = baseline_eims___15,
    ethnicity = gi_q1,
    ethnicity_detail = gi_q1_other_3,
    is_in_good_health = q_hc_goodhealth,
    not_in_good_health_detail = gi_q1_other_2,
    on_long_term_meds = q_hc_longterm_medication,
    meds_list = gi_q1_other,
    current_smoker = smokeryn1,
    smoking_at_diagnosis = q3_gi,
    recent_physical_activity = q5_gi,
    recent_physical_activity_detail = q5_gi_yes,
    drinks_alcohol = q4_gi,
    alcohol_frequency = q4_gi_yes,
    fh_ibd = family_history_pfh,
    fh_relationship = relationship_pfh,
    fh_diagnosis = fhdiagnosis_pfh_1,
    fh_other = other_family_diagnosis,
    had_tonsillectomy = mh_tonsilout1,
    tonsillectomy_date = mh_tonsil_date1,
    had_appendicectomy = mh_appendx1,
    appendicectomy_date = mh_appendix_date1,
    questionnaire_comments = comment_if_any_of_these_qu,
  )
as.matrix(colnames(cleaned_demographic_df))

cleaned_demographic_df[cleaned_demographic_df == ""] <- NA

df1 <- cleaned_demographic_df

# Clean data types now
print("Performing type conversions for demographic dataframe...")

print("Integer conversions...")
# Integer Conversions
integer_cols <- c(
  "age",
  "age_at_diagnosis"
)

df1[,integer_cols] <- lapply(df1[,integer_cols], as.integer)
print("Numeric conversions...")
# Float Conversions
df1$height <- as.numeric(df1$height)
df1$weight <- as.numeric(df1$weight)
df1$bmi <- as.numeric(df1$bmi)
df1$patient_reported_health_status <- as.numeric(df1$patient_reported_health_status)

# Date Conversions
print("Date conversions...")
date_cols <- c(
  "date_of_diagnosis",
  "data_entry_date",
  "ifx_start",
  "ifx_stop",
  "ada_start",
  "ada_stop",
  "vedo_start",
  "vedo_stop",
  "uste_start",
  "uste_stop",
  "tofa_start",
  "tofa_stop",
  "golim_start",
  "golim_stop",
  "aza_start",
  "aza_stop",
  "mp_start",
  "mp_stop",
  "mtx_start",
  "mtx_stop",
  "gi_questionnaire_date"
  )

df1[,date_cols] <- lapply(df1[,date_cols], as.Date)

print("Factor conversions...")
# Factor Conversions
df1$sex <- as.factor(df1$sex)
df1$sex <- recode(df1$sex, "1"="male", "2"="female")

df1$is_healthy_control <- as.factor(df1$is_healthy_control)
df1$is_healthy_control <- recode(df1$is_healthy_control, "0"="participant", "1"="healthy_control")

df1$group <- as.factor(df1$group)
df1$group <- recode(df1$group, "1"="CD", "2"="UC", "3"="IBDU", "4"="Non-IBD", "5"="Awaiting diagnosis", "6"="HC")

df1$patient_details_complete <- as.factor(df1$patient_details_complete)
df1$patient_details_complete <- recode(df1$patient_details_complete, "0"="incomplete", "1"="unverified", "2"="complete")

df1$ibd_background_clinician_complete <- as.factor(df1$ibd_background_clinician_complete)
df1$ibd_background_clinician_complete <- recode(df1$ibd_background_clinician_complete, "0"="incomplete", "1"="unverified", "2"="complete")

df1$gidamps_participant_questionnaire_complete <- as.factor(df1$gidamps_participant_questionnaire_complete)
df1$gidamps_participant_questionnaire_complete <- recode(df1$gidamps_participant_questionnaire_complete, "0"="incomplete", "1"="unverified", "2"="complete")

df1$baseline_recruitment_type <- as.factor(df1$baseline_recruitment_type)
df1$baseline_recruitment_type <- recode(df1$baseline_recruitment_type, "1"="endoscopy", "2"="outpatient", "3"="inpatient")

df1$montreal_extent <- as.factor(df1$montreal_extent)
df1$montreal_extent <- recode(df1$montreal_extent, "0"="E1", "1"="E2", "2"="E3")

df1$montreal_severity <- as.factor(df1$montreal_severity)
df1$montreal_severity <- recode(df1$montreal_severity, "0"="S0", "1"="S1", "2"="S2", "3"="S3")

df1$past_ibd_surgery <- as.factor(df1$past_ibd_surgery)
df1$past_ibd_surgery <- recode(df1$past_ibd_surgery, "1"="yes", "2"="no")

df1$fh_ibd <- as.factor(df1$fh_ibd)
df1$fh_ibd <- recode(df1$fh_ibd, "1"="yes", "2"="no", "3"="no")

df1$ifx <- as.factor(df1$ifx)
df1$ifx <- recode(df1$ifx, "1"="yes", "2"="no")

df1$had_appendicectomy <- as.factor(df1$had_appendicectomy)
df1$had_appendicectomy <- recode(df1$had_appendicectomy, "1"="yes", "2"="no", "0"="no")
droplevels(df1$had_appendicectomy)

df1$ethnicity <- as.factor(df1$ethnicity)
df1$ethnicity <- recode(df1$ethnicity, "1"="white_european", "2"="other")

df1$current_smoker <- as.factor(df1$current_smoker)
df1$current_smoker <- recode(df1$current_smoker, "1"="smoker", "2"="ex_smoker", "3"="non_smoker")

df1$smoking_at_diagnosis <- as.factor(df1$smoking_at_diagnosis)
df1$smoking_at_diagnosis <- recode(df1$smoking_at_diagnosis, "1"="yes", "2"="no", "3"="unsure")


factor_cols <- c(
  "is_new_diagnosis",
  "montreal_l1",
  "montreal_l2",
  "montreal_l3",
  "montreal_l4",
  "montreal_b1",
  "montreal_b2",
  "montreal_b3",
  "montreal_p",
  "eims_arthralgia",
  "eims_ank_spon",
  "eims_erythema_nodosum",
  "eims_pyoderma",
  "eims_uveitis",
  "eims_scleritis",
  "eims_sacroileitis",
  "eims_none",
  "past_medical_history",
  "ada",
  "vedo",
  "uste",
  "tofa",
  "golim",
  "aza",
  "mp",
  "mtx",
  "is_in_good_health",
  "on_long_term_meds",
  "recent_physical_activity",
  "drinks_alcohol",
  "had_tonsillectomy"
)
df1[,factor_cols] <- lapply(df1[,factor_cols], as.factor)

print("Conversions complete.")
print("Saving demographic dataframe to rds and csv...")
saveRDS(df1, file=paste0(output_dir, "demographic_df.rds"))
write.csv(df1, file=paste0(output_dir, "demographic_df.csv"), row.names = FALSE)
print("Saving complete.")
print("Demographic dataframe complete.")



print("Creating and cleaning sampling dataframe...")
# Create a clean sampling data frame
sampling_df <- df %>%
  filter(redcap_repeat_instrument == "sampling") %>%
  select(-colnames(demographic_df)[-1]) %>%
  select(-ends_with("_1000")) %>%
  select(-starts_with("hbicomplications")) %>%
  select(-starts_with("sccai_complications")) %>%
  select(-c(
    study_group,
    study_group_hc,
    email,
    consent
  ))

sampling_df[sampling_df == ""] <- NA



print("Performing type conversions...")
sampling_factor_cols <- c(
  "sampling_asa",
  "sampling_aza",
  "sampling_mp",
  "sampling_ifx",
  "sampling_ada",
  "sampling_vedo",
  "sampling_uste",
  "sampling_tofa",
  "sampling_mtx",
  "sampling_goli",
  "sampling_ciclosporin",
  "blood_experiment",
  "faecal_experiment",
  "biopsy_experiment",
  "blood_sample_collected___5",
  "blood_sample_collected___6",
  "blood_sample_collected___3",
  "blood_sample_collected___4",
  "blood_sample_collected___100",
  "blood_sample_collected___200",
  "blood_sample_collected___99",
  "blood_sample_collected___2",
  "blood_sample_collected___1",
  "blood_sample_optinal_set___1",
  "blood_sample_additional___1",
  "blood_sample_additional___2",
  "faecal_sample_collected___1",
  "faecal_sample_collected___2",
  "faecal_sample_collected___3",
  "faecal_sample_collected___4",
  "faecal_sample_collected___99",
  "biopsy_sample_collected___1",
  "biopsy_sample_collected___2",
  "biopsy_sample_collected___99",
  "patient_active_symptomyn1",
  "adalimumab_test",
  "infliximab_test",
  "vedolizumab_test",
  "ustekinumab_test",
  "steroids",
  "iv_steroids",
  "rescue_therapy",
  "antibiotics",
  "endoscopy_yn",
  "endoscopy_active",
  "radiology",
  "radiology_active",
  "endoscopy_type___1",
  "endoscopy_type___2",
  "endoscopy_type___3"
)
sampling_df[,sampling_factor_cols] <- lapply(sampling_df[,sampling_factor_cols], as.factor)

print("Renaming columns...")
cleaned_sampling_df <- sampling_df %>%
  rename(
    bloods_taken = blood_experiment,
    stools_taken = faecal_experiment,
    biopsies_taken = biopsy_experiment,
    tube_edta_10_ml = blood_sample_collected___5,
    tube_edta_genotyping = blood_sample_collected___6,
    tube_ccfdna = blood_sample_collected___3,
    tube_rna = blood_sample_collected___4,
    tube_optional_sample_set = blood_sample_collected___100,
    tube_additional_sample_set = blood_sample_collected___200,
    tube_other_sample = blood_sample_collected___99,
    tube_nhs_serum = blood_sample_collected___2,
    tube_nhs_edta = blood_sample_collected___1,
    tube_optional_citrate = blood_sample_optinal_set___1,
    tube_additional_edta = blood_sample_additional___1,
    tube_additional_cfdna = blood_sample_additional___2,
    tube_calprotectin = faecal_sample_collected___1,
    tube_qfit = faecal_sample_collected___2,
    tube_plain_stool = faecal_sample_collected___3,
    tube_omnigut = faecal_sample_collected___4,
    tube_other_stool = faecal_sample_collected___99,
    tube_biopsy_formalin = biopsy_sample_collected___1,
    tube_biopsy_rnalater = biopsy_sample_collected___2,
    tube_biopsy_other = biopsy_sample_collected___99,
    active_symptoms = patient_active_symptomyn1,
    sccai = sccai_total_calculation,
    hbi = hbi_total_calculation,
    date_of_nhs_bloods = blood_test_date_red,
    date_of_nhs_calprotectin = faecal_test_date_2,
    date_of_nhs_tdm = date_test_adali_inflixi,
    date_of_radiology = radiology_test_date,
    date_of_endoscopy = endoscopy_date,
    endoscopy_performed = endoscopy_yn,
    endoscopy_result = endoscopy_result_endcospy,
    pathology_result = histopathology_report_text,
    radiology_performed = radiology,
    ct_result = radilogy_r_ctabdo_pelvic,
    mri_sb_result = radiology_r_mri_sml_bowel,
    mri_pelvis_result = radiology_r_mri_pelvis,
    endoscopy_type_colonoscopy = endoscopy_type___2,
    endoscopy_type_flexi = endoscopy_type___3,
    endoscopy_type_other = endoscopy_type___1,
  )
  
cleaned_sampling_df <- cleaned_sampling_df %>%
  select(-starts_with("sccai_"))

print("Integer conversions...")
sampling_integer_cols <- c(
  "redcap_repeat_instance",
  "sccai",
  "hbi",
  "mayo_total",
  "uceis_total",
  "sescd_total",
  "iv_steroids_day"
)
cleaned_sampling_df[,sampling_integer_cols] <- lapply(cleaned_sampling_df[,sampling_integer_cols], as.integer)

print("Numeric conversions...")
sampling_numeric_cols <- c(
  "haemoglobin",
  "haematocrit_lab",
  "red_cell_count",
  "white_cell_count",
  "neutrophils_lab",
  "lymphocytes_lab",
  "monocytes_lab",
  "basophils_lab",
  "eosinophils_lab",
  "plt_lab",
  "urea_lab",
  "sodium_lab",
  "potassium_lab",
  "creatinine_lab",
  "egfr_lab",
  "crp",
  "albumin",
  "calprotectin",
  "drug_level_adalimumab",
  "drug_level_inflxi",
  "drug_level_vedo",
  "drug_level_uste"
)
cleaned_sampling_df[,sampling_numeric_cols] <- lapply(cleaned_sampling_df[,sampling_numeric_cols], as.numeric)

print("Date conversions...")
sampling_date_cols <- c(
  "date_of_nhs_bloods",
  "date_of_nhs_calprotectin",
  "date_of_nhs_tdm",
  "date_of_radiology",
  "date_of_endoscopy",
  "sampling_date",
  "faecal_test_date"
)
cleaned_sampling_df[,sampling_date_cols] <- lapply(cleaned_sampling_df[,sampling_date_cols], as.Date)

print("Factor conversions...")
cleaned_sampling_df$ct_abdomen_and_or_pelvis <- as.factor(cleaned_sampling_df$ct_abdomen_and_or_pelvis)
cleaned_sampling_df$ct_abdomen_and_or_pelvis <- recode(cleaned_sampling_df$ct_abdomen_and_or_pelvis, "1"="yes", "2"="no")

cleaned_sampling_df$mri_small_bowel <- as.factor(cleaned_sampling_df$mri_small_bowel)
cleaned_sampling_df$mri_small_bowel <- recode(cleaned_sampling_df$mri_small_bowel, "1"="yes", "2"="no")

cleaned_sampling_df$mri_pelvis <- as.factor(cleaned_sampling_df$mri_pelvis)
cleaned_sampling_df$mri_pelvis <- recode(cleaned_sampling_df$mri_pelvis, "1"="yes", "2"="no")

cleaned_sampling_df$ibd_status <- as.factor(cleaned_sampling_df$ibd_status)
cleaned_sampling_df$ibd_status <- recode(cleaned_sampling_df$ibd_status, "0"="Biochemical remission", "1"="Remission", "2"="Active", "3"="Highly active", "4"="Not applicable")

cleaned_sampling_df$sampling_complete <- as.factor(cleaned_sampling_df$sampling_complete)
cleaned_sampling_df$sampling_complete <- recode(cleaned_sampling_df$sampling_complete, "0"="incomplete", "1"="unverified", "2"="complete")


print("Conversions complete.")

print("Saving sampling dataframe to rds and csv...")
saveRDS(cleaned_sampling_df, file=paste0(output_dir, "sampling_df.rds"))
write.csv(cleaned_sampling_df, file=paste0(output_dir, "sampling_df.csv"), row.names = FALSE)
print("Saving complete.")
print("Sampling dataframe complete.")


print("Create merged dataset...")
merged_df <- left_join(cleaned_sampling_df, df1, by="study_id")
saveRDS(merged_df, file=paste0(output_dir, "merged_df.rds"))
write.csv(merged_df, file=paste0(output_dir, "merged_df.csv"), row.names = FALSE)
print("Saving complete.")
print("Merged dataframe complete.")

print("Script completed successfully.")

