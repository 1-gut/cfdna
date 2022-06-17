# GIDAMPS BIOANALYZER ANALYSIS
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
output_dir <- "output/gidamps/bioanalyzer/" 
input_annotation_file <- "data/gidamps/bioanalyzer_annotations/main_annotations.rds"
input_bioanalyzer_files_dir <- "data/gidamps/bioanalyzer"
# -----------------------------------------------------------------------------
print("Loading required packages...")
library("bioanalyzeR")
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggprism)
source("./theme_options.R")
print("Complete.")
# -----------------------------------------------------------------------------
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
  bioanalyzer_files[[13]]
)
print("Complete.")
# -----------------------------------------------------------------------------
# Annotate our samples with clinical data
print("Annotating data...")
annotation.file <- readRDS(file = input_annotation_file)
df <- annotate.electrophoresis(cfdna_bioanalyzer, annotation.file)
print("Complete.")
# Collapse Non-IBD into HC (deprecated)
# annotation.file$group[annotation.file$group == "Non-IBD"] <- "HC"
# df_simplified <- annotate.electrophoresis(cfdna_bioanalyzer, annotation.file)
# -----------------------------------------------------------------------------
print("Removing outliers...")
cfdna_bioanalyzer_outliers_removed <- subset(df, sample.name != 258)
cfdna_bioanalyzer_outliers_removed <- subset(cfdna_bioanalyzer_outliers_removed, sample.name != 363)
cfdna_bioanalyzer_outliers_removed <- subset(cfdna_bioanalyzer_outliers_removed, sample.name != 287)
annotation.file <- subset(annotation.file, sample.name != 258)
annotation.file <- subset(annotation.file, sample.name != 363)
annotation.file <- subset(annotation.file, sample.name != 287)
print("Complete.")
# -----------------------------------------------------------------------------
print("Begin plotting...")

# Identify outliers with this
#
# qplot.electrophoresis(
#   cfdna_bioanalyzer,
#   y = "fluorescence",
#   ylim = c(0, 450),
#   log = "x",
# )
# Outliers are 258 and 363 for the high tail fragments, 287 for the short spike

nrow(annotation.file)

count <- nrow(annotation.file)
count_label <- str_c("n=", count)

# Plot all samples
# -----------------------------------------------------------------------------
plot.all <- qplot.electrophoresis(
  cfdna_bioanalyzer_outliers_removed,
  facets = NULL,
  y = "fluorescence",
  ylim = c(0, 300),
  log = "x",
) + annotate(
  "text",
  x = 7000,
  y = 400,
  size = 8,
  label = count_label
) + 
  labs(x = "Length (bp)", y = "Fluorescence (fu)", title = "cfDNA Fragment Size (All Samples)") +
  theme_options

print(plot.all)
ggsave(paste0(output_dir, "bioanalyzer_all.png"), plot.all, dpi = 300, width = 2500, height = 2000, units = "px")

# Plot traces by group and activity
# -----------------------------------------------------------------------------
cfdna_bioanalyzer_outliers_removed_ibd_only <- subset(cfdna_bioanalyzer_outliers_removed, group == "UC" | group == "CD")
cfdna_bioanalyzer_outliers_removed_ibd_only <- subset(cfdna_bioanalyzer_outliers_removed_ibd_only, ibd_status != "")

# Retrieve counts for each plot
ibd_status_group_counts <- rename(count(annotation.file, ibd_status, group), Freq = n)
ibd_status_group_counts <- drop_na(ibd_status_group_counts)
ibd_status_group_counts$label <- str_c("n=", ibd_status_group_counts$Freq)
ibd_status_group_counts <- subset(ibd_status_group_counts, ibd_status != "Not applicable")
ibd_status_counts_without_ibdu <- subset(ibd_status_group_counts, group != "IBDU")
ibd_status_counts_without_ibdu

plot.by.group.and.activity <- qplot.electrophoresis(
  cfdna_bioanalyzer_outliers_removed_ibd_only,
  facets = group ~ ibd_status,
  y = "fluorescence",
  color = sample.name,
  show.peaks = "none",
  region.alpha = NA,
  ylim = c(0, 250),
  log = "x",
  include.ladder = FALSE,
) + 
  geom_text(
  size = 8,
  data = ibd_status_counts_without_ibdu,
  mapping = aes(x = 7000, y = 200, label = label),
  hjust = 1.05,
  vjust = 1.5
) +
  labs(x = "Length (bp)", y = "Fluorescence (fu)", title = "cfDNA Fragment Size by Activity in CD and UC") + theme_options

print(plot.by.group.and.activity)
ggsave(paste0(output_dir, "bioanalyzer_group_activity.png"), plot.by.group.and.activity, dpi = 300, width = 5000, height = 2800, units = "px")

# Plot traces by activity
# -----------------------------------------------------------------------------

# Create n=x count labels for annotating the plot
annotation.uc_cd <- subset(annotation.file, group == "CD" | group == "UC")
ibd_status_counts <- as.data.frame(table(annotation.uc_cd$ibd_status))
ibd_status_counts <- rename(ibd_status_counts, ibd_status = Var1, label = Freq)
ibd_status_counts$label <- as.character(ibd_status_counts$label)
ibd_status_counts$label <- str_c("n=", ibd_status_counts$label)
ibd_status_counts <- subset(ibd_status_counts, ibd_status != "Not applicable")
print(ibd_status_counts)

# Plot by activity
plot.by.ibd.activity <- qplot.electrophoresis(
  cfdna_bioanalyzer_outliers_removed_ibd_only,
  facets = . ~ ibd_status,
  y = "fluorescence",
  color = sample.name,
  show.peaks = "none",
  region.alpha = NA,
  ylim = c(0, 250),
  log = "x",
  include.ladder = FALSE,
) + 
  geom_text(
  size = 8,
  data = ibd_status_counts,
  mapping = aes(x = 7000, y = 200, label = label),
  hjust = 1.05,
  vjust = 1.5
) + 
  labs(x = "Length (bp)", y = "Fluorescence (fu)", title = "cfDNA Fragment Size by IBD Activity") + 
  theme_options

print(plot.by.ibd.activity) 
ggsave(paste0(output_dir, "bioanalyzer_activity.png"), plot.by.ibd.activity, dpi = 300, width = 5000, height = 1500, units = "px")

# Plot traces by group
# -----------------------------------------------------------------------------
# Retrieve the 4 main groups
group_df <- subset(cfdna_bioanalyzer_outliers_removed, group == "CD" | group == "UC" | group == "NC" | group == "HC")

# Create n=x count labels for annotating the plot
group_counts <- as.data.frame(table(annotation.file$group))
group_counts <- rename(group_counts, group = Var1, label = Freq)
group_counts$label <- as.character(group_counts$label)
group_counts$label <- str_c("n=", group_counts$label)
group_counts <- subset(group_counts, group == "CD" | group == "UC" | group == "NC" | group == "HC")
print(group_counts)

# Plot by group
plot.by.group <- qplot.electrophoresis(
  group_df,
  facets = . ~ group,
  y = "fluorescence",
  color = sample.name,
  show.peaks = "none",
  region.alpha = NA,
  ylim = c(0, 250),
  log = "x",
  include.ladder = FALSE,
) + geom_text(
  size = 8,
  data = group_counts,
  mapping = aes(x = 7000, y = 200, label = label),
  hjust = 1.05,
  vjust = 1.5
) +
  labs(x = "Length (bp)", y = "Fluorescence (fu)", title = "cfDNA Fragment Size by Group") +
  theme_options

# Print & Save
print(plot.by.group)
ggsave(paste0(output_dir, "bioanalyzer_group.png"), plot.by.group, dpi = 300, width = 5000, height = 1500, units = "px")
# -----------------------------------------------------------------------------
print("Script completed successfully.")
