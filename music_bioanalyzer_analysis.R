# GIDAMPS BIOANALYZER ANALYSIS
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
output_dir <- "output/music/bioanalyzer/" 
input_annotation_file <- "data/music/bioanalyzer/music_bioanalyzer_annotations_manual.csv"
input_music_clinical_data <- "data/music/clinical/music_clinical_2022-07-15.rds"
input_mito_trajectory_file <- "data/music/bioanalyzer/music_rising_mito_trajectory.csv"
input_endoscopic_healing_file <- "data/music/bioanalyzer/music_mucosal_healing.csv"
input_bioanalyzer_files_dir <- "data/music/bioanalyzer"
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
  bioanalyzer_files[[10]]
)
print("Complete.")
# -----------------------------------------------------------------------------
# Annotate our samples with clinical data
print("Annotating data...")
annotation.file <- read.csv(file = input_annotation_file)
# DEPENDENCY: Requires input MUSIC clinical data from music_get_clinical_data.R
# Run the above file first.
clinical_data <- readRDS(input_music_clinical_data)

cleaned_annotation_file <- dplyr::left_join(annotation.file, clinical_data, by=c("study_id", "redcap_event_name"))
levels(cleaned_annotation_file$study_group_name) <- c(levels(cleaned_annotation_file$study_group_name), "NC")
cleaned_annotation_file$study_group_name[is.na(cleaned_annotation_file$study_group_name)] <- "NC"
cleaned_annotation_file$redcap_event_name <- as.factor(cleaned_annotation_file$redcap_event_name)

# Add mito trajectory data here
rising_mito_csv <- read.csv(input_mito_trajectory_file)
cleaned_annotation_file <- dplyr::left_join(cleaned_annotation_file, rising_mito_csv, by="study_id")
cleaned_annotation_file$rising_mitochondrial_cfdna[is.na(cleaned_annotation_file$rising_mitochondrial_cfdna)] <- "Falling trajectory"
cleaned_annotation_file$rising_mitochondrial_cfdna[cleaned_annotation_file$study_group_name=="NC"] <- NA

# Add endoscopic healing data here (manually generated)
endoscopic_healing <- read.csv(input_endoscopic_healing_file)
cleaned_annotation_file <- dplyr::left_join(cleaned_annotation_file, endoscopic_healing, by="study_id")
cleaned_annotation_file$complete_mucosal_healing[cleaned_annotation_file$complete_mucosal_healing==""] <- NA
cleaned_annotation_file$endoscopic_improvement[cleaned_annotation_file$endoscopic_improvement==""] <- NA

# Rename redcap event name timepoints
cleaned_annotation_file$redcap_event_name <- recode_factor(cleaned_annotation_file$redcap_event_name,
                                                           timepoint_1 = "Baseline",
                                                           timepoint_2 = "3 months",
                                                           timepoint_3 = "6 months")

df_interim <- annotate.electrophoresis(cfdna_bioanalyzer, cleaned_annotation_file)
df <- subset(df_interim, to_discard != "yes")
cleaned_annotation_file <- cleaned_annotation_file %>% subset(to_discard != "yes")


print("Complete.")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
print("Begin plotting...")

# See all plots
# qplot.electrophoresis(df, y = "fluorescence", log = "x")
# nrow(annotation.file)
# 
count <- nrow(cleaned_annotation_file)
count_label <- str_c("n=", count)

# Plot all samples
# -----------------------------------------------------------------------------
plot_all <- qplot.electrophoresis(
  df,
  facets = NULL,
  y = "fluorescence",
  ylim = c(0, 200),
  log = "x",
) + annotate(
  "text",
  x = 7000,
  y = 150,
  size = 8,
  label = count_label
) + 
  labs(x = "Length (bp)", y = "Fluorescence (fu)", title = "cfDNA Fragment Size (MUSIC)") +
  theme_options

print(plot_all)
ggsave(paste0(output_dir, "bioanalyzer_all.png"), plot_all, dpi = 300, width = 2500, height = 2000, units = "px")


# Plot traces by group
# -----------------------------------------------------------------------------

# Create n=x count labels for annotating the plot
group_counts <- rename(count(cleaned_annotation_file, study_group_name), Freq = n)
group_counts$label <- str_c("n=", group_counts$Freq)

# Plot by group
plot_by_group <- qplot.electrophoresis(
  df,
  facets = . ~ study_group_name,
  y = "fluorescence",
  color = sample.name,
  show.peaks = "none",
  region.alpha = NA,
  ylim = c(0, 150),
  log = "x",
  include.ladder = FALSE,
) + geom_text(
  size = 8,
  data = group_counts,
  mapping = aes(x = 7000, y = 100, label = label),
  hjust = 1.05,
  vjust = 1.5
) +
  labs(x = "Length (bp)", y = "Fluorescence (fu)", title = "cfDNA Fragment Size by Group (MUSIC)") +
  theme_options

# Print & Save
print(plot_by_group)
ggsave(paste0(output_dir, "bioanalyzer_by_group.png"), plot_by_group, dpi = 300, width = 5000, height = 1500, units = "px")
# -----------------------------------------------------------------------------




# Plot traces by group and activity
# -----------------------------------------------------------------------------
df_without_nc <- subset(df, study_group_name != "NC")
annotations_without_nc <- subset(cleaned_annotation_file, study_group_name != "NC")

# Retrieve counts for each plot
group_and_time_counts <- rename(count(annotations_without_nc, study_group_name, redcap_event_name), Freq = n)
group_and_time_counts$label <- str_c("n=", group_and_time_counts$Freq)


plot_by_group_and_time <- qplot.electrophoresis(
  df_without_nc,
  facets = study_group_name ~ redcap_event_name,
  y = "fluorescence",
  color = sample.name,
  show.peaks = "none",
  region.alpha = NA,
  ylim = c(0, 150),
  log = "x",
  include.ladder = FALSE,
) + 
  geom_text(
  size = 8,
  data = group_and_time_counts,
  mapping = aes(x = 7000, y = 100, label = label),
  hjust = 1.05,
  vjust = 1.5
) +
  labs(x = "Length (bp)", y = "Fluorescence (fu)", title = "cfDNA Fragment Size by Group Over Time (MUSIC)") + theme_options

print(plot_by_group_and_time)
ggsave(paste0(output_dir, "bioanalyzer_by_group_over_time.png"), plot_by_group_and_time, dpi = 300, width = 5000, height = 2800, units = "px")
# -----------------------------------------------------------------------------


# Plot traces by mitochondrial trajectory
# -----------------------------------------------------------------------------

# Create n=x count labels for annotating the plot
mito_counts <- rename(count(annotations_without_nc, rising_mitochondrial_cfdna, redcap_event_name), Freq = n)
mito_counts$label <- str_c("n=", mito_counts$Freq)

# Plot by group
plot_by_mito_trajectory <- qplot.electrophoresis(
  df_without_nc,
  facets = rising_mitochondrial_cfdna ~ redcap_event_name,
  y = "fluorescence",
  color = sample.name,
  show.peaks = "none",
  region.alpha = NA,
  ylim = c(0, 150),
  log = "x",
  include.ladder = FALSE,
) + geom_text(
  size = 8,
  data = mito_counts,
  mapping = aes(x = 7000, y = 100, label = label),
  hjust = 1.05,
  vjust = 1.5
) +
  labs(x = "Length (bp)", y = "Fluorescence (fu)", title = "cfDNA Fragment Size by Mitochondrial cfDNA Trajectory") +
  theme_options

# Print & Save
print(plot_by_mito_trajectory)
ggsave(paste0(output_dir, "bioanalyzer_by_mito_trajectory.png"), plot_by_mito_trajectory, dpi = 300, width = 5000, height = 2800, units = "px")
# -----------------------------------------------------------------------------

# Plot traces by complete mucosal healing
# -----------------------------------------------------------------------------
df_mucosal_healing <- subset(df_without_nc, !is.na(complete_mucosal_healing))
annotations_mucosal_healing <- subset(annotations_without_nc, !is.na(complete_mucosal_healing))


# Create n=x count labels for annotating the plot
mucosal_healing_counts <- rename(count(annotations_mucosal_healing, complete_mucosal_healing, redcap_event_name), Freq = n)
mucosal_healing_counts$label <- str_c("n=", mucosal_healing_counts$Freq)

# Plot by group
plot_by_mucosal_healing <- qplot.electrophoresis(
  df_mucosal_healing,
  facets = complete_mucosal_healing ~ redcap_event_name,
  y = "fluorescence",
  color = sample.name,
  show.peaks = "none",
  region.alpha = NA,
  ylim = c(0, 150),
  log = "x",
  include.ladder = FALSE,
) + geom_text(
  size = 8,
  data = mucosal_healing_counts,
  mapping = aes(x = 7000, y = 100, label = label),
  hjust = 1.05,
  vjust = 1.5
) +
  labs(x = "Length (bp)", y = "Fluorescence (fu)", title = "cfDNA Fragment Size by Complete Mucosal Healing") +
  theme_options

# Print & Save
print(plot_by_mucosal_healing)
ggsave(paste0(output_dir, "bioanalyzer_by_mucosal_healing.png"), plot_by_mucosal_healing, dpi = 300, width = 5000, height = 2800, units = "px")
# -----------------------------------------------------------------------------

# Plot traces by endoscopic improvement
# -----------------------------------------------------------------------------
df_endoscopic_improvement <- subset(df_without_nc, !is.na(endoscopic_improvement))
annotations_endoscopic_improvement <- subset(annotations_without_nc, !is.na(endoscopic_improvement))


# Create n=x count labels for annotating the plot
endoscopic_improvement_counts <- rename(count(annotations_endoscopic_improvement, endoscopic_improvement, redcap_event_name), Freq = n)
endoscopic_improvement_counts$label <- str_c("n=", endoscopic_improvement_counts$Freq)

# Plot by group
plot_by_endoscopic_improvement <- qplot.electrophoresis(
  df_endoscopic_improvement,
  facets = endoscopic_improvement ~ redcap_event_name,
  y = "fluorescence",
  color = sample.name,
  show.peaks = "none",
  region.alpha = NA,
  ylim = c(0, 150),
  log = "x",
  include.ladder = FALSE,
) + geom_text(
  size = 8,
  data = endoscopic_improvement_counts,
  mapping = aes(x = 7000, y = 100, label = label),
  hjust = 1.05,
  vjust = 1.5
) +
  labs(x = "Length (bp)", y = "Fluorescence (fu)", title = "cfDNA Fragment Size by Endoscopic Improvement") +
  theme_options

# Print & Save
print(plot_by_endoscopic_improvement)
ggsave(paste0(output_dir, "bioanalyzer_by_endoscopic_improvement.png"), plot_by_endoscopic_improvement, dpi = 300, width = 5000, height = 2800, units = "px")
# -----------------------------------------------------------------------------



print("Script completed successfully.")
