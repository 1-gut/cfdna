# Runs all analysis in this directory
# R version 4.2.0
print("Retrieving clinical data from redcap...")
source("./gidamps_get_clinical_data.R")
source("./music_get_clinical_data.R")
print("Data retrieval completed.")

print("Run dPCR and demographics...")
source("./gidamps_dpcr_analysis.R")


print("Script completed successfully.")