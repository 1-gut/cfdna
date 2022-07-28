# GIDAMPS CLINICAL ANALYSIS
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
output_dir <- "output/gidamps/clinical/"
input_gidamps_merged_df <- "data/gidamps/clinical/merged_df.rds"
# -----------------------------------------------------------------------------
# Description
# -----------------------------------------------------------------------------
# Loads dPCR data from prepared excel file and merges with clinical data from
# music_get_clinical_data.R
# Creates graphs into output directory.
# -----------------------------------------------------------------------------
print("Loading required packages...")
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(summarytools)
source("./theme_options.R")
print("Complete.")

df <- readRDS(input_gidamps_merged_df)


add_center_name <- function(df){
  study_id = df[1]
  if(substr(study_id, 1,4) == "136-") {
    center_name <- "Glasgow"
  } else if (substr(study_id,1,4) == "138-") {
    center_name <- "Dundee"
  } else {
    center_name <- "Edinburgh"
  }
  return(center_name)
}

df$center_name <- apply(df, 1, add_center_name)

glasgow <- df %>% subset(center_name == "Glasgow")
edinburgh <- df %>% subset(center_name == "Edinburgh")
dundee <- df %>% subset(center_name == "Dundee")

view(dfSummary(glasgow), file=paste0(output_dir, "glasgow.html"))
view(dfSummary(edinburgh), file=paste0(output_dir, "edinburgh.html"))
view(dfSummary(dundee), file=paste0(output_dir, "dundee.html"))
