# GIDAMPS MAKE BIOANALYZER ANNOTATIONS
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
output_dir <- "data/gidamps/bioanalyzer_annotations/"
gidamps_merged_clinical_file <- "data/gidamps/clinical/merged_df.rds" # rds file
sample_id_date_list <- "data/gidamps/bioanalyzer_annotations/sample_date_15032022.csv"

input_bioanalyzer_files_dir <- "data/gidamps/bioanalyzer"
# -----------------------------------------------------------------------------
# Description
# -----------------------------------------------------------------------------
# Loads dPCR data from prepared excel file and merges with clinical data from
# music_get_clinical_data.R
# Creates graphs into output directory.
# -----------------------------------------------------------------------------
print("Loading required packages...")
library("bioanalyzeR")
library(dplyr)
library(stringr)
library(ggplot2)
print("Complete.")

# Load data files from data directory (need to manually import each file)
print("Loading and reading data...")
bioanalyzer_files <- list.files(input_bioanalyzer_files_dir, full.names = TRUE)
cfdna_bioanalyzer <- read.electrophoresis(
  bioanalyzer_files[[1]],
  bioanalyzer_files[[2]],
  bioanalyzer_files[[3]],
  bioanalyzer_files[[4]],
  bioanalyzer_files[[5]],
  bioanalyzer_files[[6]],
  bioanalyzer_files[[7]],
  bioanalyzer_files[[8]],
  bioanalyzer_files[[9]],
  bioanalyzer_files[[10]],
  bioanalyzer_files[[11]],
  bioanalyzer_files[[12]],
  bioanalyzer_files[[13]],
  bioanalyzer_files[[14]]
)
print("Complete.")

# Extract the study_ids 
print("Extracting study IDs from bioanalyzer files...")
sample_list_df <- cfdna_bioanalyzer[["samples"]]
df = subset(sample_list_df, select = c("sample.name") )
df2 = subset(df, sample.name!="Ladder")
df2$study_id <- as.numeric(str_extract(df2$sample.name, "[0-9]+"))

# Load sample_date data from annotations (this came from original cfDNA data file)
print("Loading sample_date from original cfDNA extraction excel...")
sample_dates <- read.csv(sample_id_date_list, fileEncoding = "UTF-8-BOM")
sample_dates$study_id <- as.numeric(str_extract(sample_dates$Sample.number, "[0-9]+"))
sample_dates$Sample.date <- as.Date(sample_dates$Sample.date, format="%d/%m/%Y")
sample_dates <- subset(sample_dates, select=c("study_id", "Sample.date"))

# Join sample date onto df2 and rename it nicely for merging
print("Join sample dates onto bioanalyzer samples...")
df3 = left_join(df2, sample_dates, by="study_id")
df3 = rename(df3, sampling_date = Sample.date)

# Read GID dataframe in and format it for merging
print("Read GI-DAMPs data in...")
gid <- readRDS(file=gidamps_merged_clinical_file)
gid$study_id <- as.numeric(gid$study_id) # BEWARE! This step drops all the Glasgow patients 136-etc.
gid$sampling_date <- as.Date(gid$sampling_date, format="%Y-%m-%d")

print("Creating annotations file...")
annotations_df <- left_join(df3, gid, by=c("study_id", "sampling_date"))
# Clean up annotations df
print("Cleaning up annotations file...")
dedup <- annotations_df[!duplicated(annotations_df), ]

# Label the negative controls
levels(dedup$group) <- c(levels(dedup$group), "NC")
dedup$group[dedup$study_id>10000 & is.na(dedup$group)] <- "NC"

print("Saving annotations file...")
saveRDS(dedup, file=paste0(output_dir, "main_annotations.rds"))
write.csv(dedup, file=paste0(output_dir, "main_annotations.csv"))
print("Script completed.")
