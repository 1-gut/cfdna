# DPCR ASSAY REPRODUCIBILITY STATS
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
output_dir <- "output/dpcr_assay/" 
input_positive_control_repeated_measurements <- "data/dpcr_assay/cfdna_repeatability_08062022.csv"
# -----------------------------------------------------------------------------
library("MBESS")
df <- read.csv(input_positive_control_repeated_measurements)

sink(file=paste(output_dir, "dpcr_assay_stats.txt"))
cat("COX3 REPRODUCIBILITY STATS\n")
MBESS::ci.cv(data=df$cox3_repeatability, conf.level=.95)
cat("ND2 REPRODUCIBILITY STATS\n")
MBESS::ci.cv(data=df$nd2_repeatability, conf.level=.95)
sink()
print("Script completed.")
