# MUSIC DPCR ANALYSIS
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------

# RDS file
input_music_clinical_data <- "data/music/clinical/music_clinical_2022-06-10.rds"

input_mito_trajectory_file <- "data/music/music_rising_mito_trajectory.csv"
input_endoscopic_healing_file <- "data/music/music_mucosal_healing.csv"

# Reads from sheet name 'cleaned_data'
input_dpcr_file <- "data/music/dpcr/music_cfdna_master_02062022.xlsx"
output_dir <- "output/music/dpcr/"
cohort_dir <- "output/music/cohort/"
# -----------------------------------------------------------------------------
# Description
# -----------------------------------------------------------------------------
# Loads dPCR data from prepared excel file and merges with clinical data from
# music_get_clinical_data.R
# Creates graphs into output directory.
# -----------------------------------------------------------------------------
print("Loading required packages...")
library(dplyr)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(summarytools)
library(compareGroups)
library(forcats)
source("./theme_options.R")
print("Complete.")

# Loading Data
print("Loading data...")
dpcr <- read_excel(input_dpcr_file, sheet = "cleaned_data")
music <- readRDS(input_music_clinical_data)
print("Complete.")

# Merge GIDAMPs onto dPCR Data
print("Creating merged file...")
df <- left_join(dpcr, music, by = c("study_id", "redcap_event_name"))
print("Complete.")

# Add endoscopic healing and mito trajectory groups
df_mito <- read.csv(input_mito_trajectory_file)
df <- left_join(df, df_mito, by="study_id")
df$rising_mitochondrial_cfdna[is.na(df$rising_mitochondrial_cfdna)] <- "Non-rising trajectory"
df$rising_mitochondrial_cfdna <- factor(df$rising_mitochondrial_cfdna, levels=c("Rising trajectory", "Non-rising trajectory"))

df_endoscopy <- read.csv(input_endoscopic_healing_file)
df <- left_join(df, df_endoscopy, by="study_id")
df$complete_mucosal_healing[df$complete_mucosal_healing==""] <- NA
df$endoscopic_improvement[df$endoscopic_improvement==""] <- NA
df$complete_mucosal_healing <- factor(df$complete_mucosal_healing, levels=c("Yes", "No"))
df$endoscopic_improvement <- factor(df$endoscopic_improvement, levels=c("Yes", "No"))
# Create simplified df for analysis
# Choose your variables here.

df$cox3_log <- log10(df$cox3_sapphire)
df$nd2_log <- log10(df$nd2_sapphire)

# Fix UCEIS and SESCD

df$uceis_final <- ifelse(
  is.na(df$uceis_total_score),
  df$uceis_noncalc,
  df$uceis_total_score
)
uceis_check <- df[, c("uceis_final", "uceis_total_score", "uceis_noncalc")]

df$sescd_final <- ifelse(
  is.na(df$cdeis_total_score),
  df$sescd_noncalc,
  df$cdeis_total_score
)
sescd_check <- df[, c("sescd_final", "cdeis_total_score", "sescd_noncalc")]

timepoint_labels <- c(
  "timepoint_1" = "Baseline",
  "timepoint_2" = "3 months",
  "timepoint_3" = "6 months",
  "timepoint_4" = "9 months",
  "timepoint_5" = "12 months"
)

give.n <- function(x) {
  return(c(y = median(x) + 0.3, label = length(x)))
}


# ====================================================================
# Cohort Analysis
# ====================================================================

df_baseline <- subset(df, redcap_event_name == "timepoint_1")
df_baseline_cd <- subset(df_baseline, study_group_name == "CD")
df_baseline_uc <- subset(df_baseline, study_group_name == "UC")

view(dfSummary(df_baseline),
  file = paste0(cohort_dir, "music_dpcr_cohort_all.html")
)
view(dfSummary(df_baseline_cd),
  file = paste0(cohort_dir, "music_dpcr_cohort_cd.html")
)
view(dfSummary(df_baseline_uc),
  file = paste0(cohort_dir, "music_dpcr_cohort_uc.html")
)


sink(file = paste0(cohort_dir, "cohort_statistics.txt"))
demographics <- compareGroups(study_group_name ~
  age +
  sex +
  bmi_weight +
  center +
  haemoglobin +
  red_cell_count +
  white_cell_count +
  plt_lab +
  albumin +
  crp +
  calprotectin +
  uceis_final +
  sescd_final +
  cucq_total +
  physician_global_assessment +
  sccai_total +
  mayo_total +
  hbi_total,
data = df_baseline
)
createTable(demographics)
sink()
# ====================================================================
# Longitudinal Boxplots
# ====================================================================

plot.longitudinal.boxplot <- function(y, y_label) {

  # Plots boxplots of chosen variable for all patients, then UC then Crohn's.

  file_string <- paste0("boxplot_", y, ".png")
  file_string.uc <- paste0("boxplot_", y, "_uc.png")
  file_string.cd <- paste0("boxplot_", y, "_cd.png")
  file_string.combined <- paste0("boxplot_", y, "_combined.png")

  p <- df %>%
    ggplot(aes(
      x = redcap_event_name, y = .data[[y]],
      fill = redcap_event_name
    )) +
    stat_boxplot(geom = "errorbar", width = 0.3) +
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.05) +
    xlab("Timepoints") +
    ylab(y_label) +
    labs(
      title = paste0(y_label, " over time (all patients)")
    ) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_x_discrete(labels = timepoint_labels) +
    theme_options
  print(p)
  ggplot2::ggsave(paste0(output_dir, file_string), p,
    dpi = 300, width = 2000, height = 1500, units = "px"
  )

  p.uc <- df %>%
    subset(study_group_name == "UC") %>%
    ggplot(aes(
      x = redcap_event_name, y = .data[[y]],
      fill = redcap_event_name
    )) +
    stat_boxplot(geom = "errorbar", width = 0.3) +
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.05) +
    xlab("Timepoints") +
    ylab(y_label) +
    labs(
      title = "Ulcerative Colitis"
    ) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_x_discrete(labels = timepoint_labels) +
    theme_options
  print(p.uc)
  ggplot2::ggsave(paste0(output_dir, file_string.uc), p.uc,
    dpi = 300, width = 2000, height = 1500, units = "px"
  )

  p.cd <- df %>%
    subset(study_group_name == "CD") %>%
    ggplot(aes(
      x = redcap_event_name, y = .data[[y]],
      fill = redcap_event_name
    )) +
    stat_boxplot(geom = "errorbar", width = 0.3) +
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.05) +
    xlab("Timepoints") +
    ylab(y_label) +
    labs(
      title = "Crohn's Disease"
    ) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_x_discrete(labels = timepoint_labels) +
    theme_options
  print(p.cd)
  ggplot2::ggsave(paste0(output_dir, file_string.cd), p.cd,
    dpi = 300, width = 2000, height = 1500, units = "px"
  )

  p.combined <- gridExtra::grid.arrange(p.cd, p.uc, nrow = 1)
  ggplot2::ggsave(paste0(output_dir, file_string.combined), p.combined,
    dpi = 300, width = 4000, height = 1500, units = "px"
  )
}

# plot.longitudinal.boxplot(y="cox3_log", y_label="COX3")
# plot.longitudinal.boxplot(y="nd2_log", y_label="ND2")
# plot.longitudinal.boxplot(y="total_cfdna", y_label="Total cfDNA")
plot.longitudinal.boxplot(y="crp", y_label="CRP")
plot.longitudinal.boxplot(y="calprotectin", y_label="Calprotectin")
plot.longitudinal.boxplot(y="cucq_total", y_label="CUCQ-32")
plot.longitudinal.boxplot(y="hbi_total", y_label="HBI")
plot.longitudinal.boxplot(y="sccai_total", y_label="SCCAI")
plot.longitudinal.boxplot(y="mayo_total", y_label="Mayo Clinical Score")


# Special UC Clinical Score Plot
p.hbi <- df %>%
  subset(study_group_name == "CD") %>%
  ggplot(aes(
    x = redcap_event_name, y = hbi_total,
    fill = redcap_event_name
  )) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.05) +
  xlab("Timepoints") +
  ylab("Harvey-Bradshaw Index") +
  labs(
    title = "HBI (Crohn's Disease)"
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = timepoint_labels) +
  theme_options
print(p.hbi)
ggplot2::ggsave(paste0(output_dir, "hbi_clinical_scores.png"), p.hbi,
                dpi = 300, width = 2000, height = 1500, units = "px"
)


p.sccai <- df %>%
  subset(study_group_name == "UC") %>%
  ggplot(aes(
    x = redcap_event_name, y = sccai_total,
    fill = redcap_event_name
  )) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.05) +
  xlab("Timepoints") +
  ylab("SCCAI") +
  labs(
    title = "SCCAI (Ulcerative Colitis)"
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = timepoint_labels) +
  theme_options
print(p.sccai)

p.mayo <- df %>%
  subset(study_group_name == "UC") %>%
  ggplot(aes(
    x = redcap_event_name, y = mayo_total,
    fill = redcap_event_name
  )) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.05) +
  xlab("Timepoints") +
  ylab("Mayo Clinical Score") +
  labs(
    title = "Mayo Clinical Score (Ulcerative Colitis)"
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = timepoint_labels) +
  theme_options
print(p.mayo)


p.uc_clinical <- gridExtra::grid.arrange(p.sccai, p.mayo, nrow = 1)
ggplot2::ggsave(paste0(output_dir, "uc_clinical_scores.png"), p.uc_clinical,
                dpi = 300, width = 4000, height = 1500, units = "px"
)



# ====================================================================
# Trajectory Plots
# ====================================================================

plot.trajectory <- function(y, y_label) {
  # Plots trajectory of chosen variable for all patients, then UC then Crohn's.

  file_string <- paste0("trajectory_", y, ".png")
  file_string.uc <- paste0("trajectory_", y, "_uc.png")
  file_string.cd <- paste0("trajectory_", y, "_cd.png")
  file_string.combined <- paste0("trajectory_", y, "_combined.png")

  p <- df %>%
    ggplot(aes(
      x = redcap_event_name, y = .data[[y]],
      group = study_id, col = study_id
    )) +
    geom_point(size = 3) +
    geom_line() +
    geom_label(aes(label = study_id), data = df %>%
      filter(redcap_event_name == "timepoint_3")) +
    xlab("Timepoints") +
    ylab(y_label) +
    labs(
      title = paste0(y_label, " over time (all patients)")
    ) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_x_discrete(labels = timepoint_labels) +
    theme_options
  print(p)
  ggplot2::ggsave(paste0(output_dir, file_string), p,
    dpi = 300, width = 2000, height = 1500, units = "px"
  )

  p.uc <- df %>%
    subset(study_group_name == "UC") %>%
    ggplot(aes(
      x = redcap_event_name,
      y = .data[[y]],
      group = study_id,
      col = study_id
    )) +
    geom_point(size = 3) +
    geom_line() +
    geom_label(aes(label = study_id),
      data = df %>%
        subset(study_group_name == "UC") %>%
        filter(redcap_event_name == "timepoint_3")
    ) +
    xlab("Timepoints") +
    ylab(y_label) +
    labs(
      title = "Ulcerative Colitis"
    ) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_x_discrete(labels = timepoint_labels) +
    theme_options
  # print(p.uc)
  # ggplot2::ggsave(paste0(output_dir, file_string.uc), p.uc,
  # dpi = 300, width = 2500, height = 2000, units = "px")

  p.cd <- df %>%
    subset(study_group_name == "CD") %>%
    ggplot(aes(
      x = redcap_event_name, y = .data[[y]],
      group = study_id, col = study_id
    )) +
    geom_point(size = 3) +
    geom_line() +
    geom_label(aes(label = study_id),
      data = df %>%
        subset(study_group_name == "CD") %>%
        filter(redcap_event_name == "timepoint_3")
    ) +
    xlab("Timepoints") +
    ylab(y_label) +
    labs(
      title = "Crohn's Disease"
    ) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_x_discrete(labels = timepoint_labels) +
    theme_options
  # print(p.cd)
  # ggplot2::ggsave(paste0(output_dir, file_string.cd), p.cd,
  # dpi = 300, width = 2500, height = 2000, units = "px")

  p_combined <- gridExtra::grid.arrange(p.cd, p.uc, nrow = 1)
  ggplot2::ggsave(paste0(output_dir, file_string.combined), p_combined,
    dpi = 300, width = 4000, height = 1500, units = "px"
  )
}

plot.trajectory(y = "cox3_sapphire", y_label = "COX3")
# plot.trajectory(y="cox3_log", y_label="Log COX3")
plot.trajectory(y = "nd2_sapphire", y_label = "ND2")
# plot.trajectory(y="nd2_log", y_label="Log ND2")
plot.trajectory(y = "total_cfdna", y_label = "Total cfDNA")
plot.trajectory(y = "crp", y_label = "CRP")
plot.trajectory(y = "calprotectin", y_label = "Calprotectin")
plot.trajectory(y = "cucq_total", y_label = "CUCQ-32")

# Separate COX3 against trajectory

p <- df %>%
  ggplot(aes(
    x = redcap_event_name, y = cox3_log,
    group = study_id, col = study_id
  )) +
  geom_point(size = 3) +
  geom_line() +
  geom_label(aes(label = study_id), data = df %>%
               filter(redcap_event_name == "timepoint_3")) +
  xlab("Timepoints") +
  ylab("COX3") +
  labs(title = "COX3 (Log) Trajectories") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = timepoint_labels) +
  facet_grid(
    . ~ rising_mitochondrial_cfdna
  ) +
  theme_options
print(p)
ggplot2::ggsave(paste0(output_dir, "cox3_separated_by_trajectory.png"), p,
                dpi = 300, width = 4000, height = 1500, units = "px"
)

p <- df %>%
  ggplot(aes(
    x = redcap_event_name, y = nd2_log,
    group = study_id, col = study_id
  )) +
  geom_point(size = 3) +
  geom_line() +
  geom_label(aes(label = study_id), data = df %>%
               filter(redcap_event_name == "timepoint_3")) +
  xlab("Timepoints") +
  ylab("ND2") +
  labs(title = "ND2 (Log) Trajectories") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = timepoint_labels) +
  facet_grid(
    . ~ rising_mitochondrial_cfdna
  ) +
  theme_options
print(p)
ggplot2::ggsave(paste0(output_dir, "nd2_separated_by_trajectory.png"), p,
                dpi = 300, width = 4000, height = 1500, units = "px"
)

# Plot COX3 Log against mucosal healing outcomes

p <- df %>%
  filter(!is.na(complete_mucosal_healing)) %>%
  ggplot(aes(
    x = redcap_event_name, y = cox3_log,
    group = study_id, col = study_id
  )) +
  geom_point(size = 3) +
  geom_line() +
  geom_label(aes(label = study_id), data = df %>% filter(!is.na(complete_mucosal_healing)) %>%
               filter(redcap_event_name == "timepoint_3")) +
  xlab("Timepoints") +
  ylab("COX3 (Log)") +
  labs(title = "COX3 by Complete Mucosal Healing") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = timepoint_labels) +
  facet_grid(
    . ~ complete_mucosal_healing
  ) +
  theme_options
print(p)
ggplot2::ggsave(paste0(output_dir, "cox3_by_mucosal_healing.png"), p,
                dpi = 300, width = 5000, height = 2000, units = "px"
)


p <- df %>%
  filter(!is.na(complete_mucosal_healing)) %>%
  ggplot(aes(
    x = redcap_event_name, y = nd2_sapphire,
    group = study_id, col = study_id
  )) +
  geom_point(size = 3) +
  geom_line() +
  geom_label(aes(label = study_id), data = df %>% filter(!is.na(complete_mucosal_healing)) %>%
               filter(redcap_event_name == "timepoint_3")) +
  xlab("Timepoints") +
  ylab("ND2 (Log)") +
  labs(title = "ND2 by Complete Mucosal Healing") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = timepoint_labels) +
  facet_grid(
    . ~ complete_mucosal_healing
  ) +
  theme_options
print(p)
ggplot2::ggsave(paste0(output_dir, "nd2_by_mucosal_healing.png"), p,
                dpi = 300, width = 4000, height = 1500, units = "px"
)


p <- df %>%
  filter(!is.na(complete_mucosal_healing)) %>%
  ggplot(aes(
    x = redcap_event_name, y = total_cfdna,
    group = study_id, col = study_id
  )) +
  geom_point(size = 3) +
  geom_line() +
  geom_label(aes(label = study_id), data = df %>% filter(!is.na(complete_mucosal_healing)) %>%
               filter(redcap_event_name == "timepoint_3")) +
  xlab("Timepoints") +
  ylab("Total cfDNA (ng/uL)") +
  labs(title = "Total cfDNA by Complete Mucosal Healing") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = timepoint_labels) +
  facet_grid(
    . ~ complete_mucosal_healing
  ) +
  theme_options
print(p)
ggplot2::ggsave(paste0(output_dir, "total_cfdna_by_mucosal_healing.png"), p,
                dpi = 300, width = 4000, height = 1500, units = "px"
)


# ====================================================================
# COX3 log by Gender, PGA
# ====================================================================

df %>%
  ggplot(aes(x = sex, y = cox3_log, fill = sex)) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.05) +
  xlab("Gender") +
  ylab("COX3 (copies/uL)") +
  labs(
    title = "COX3 by Gender (all patients)"
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_options
ggplot2::ggsave(paste0(output_dir, "cox3_vs_gender.png"),
  dpi = 300, width = 2000, height = 1500, units = "px"
)

df %>%
  ggplot(aes(x = physician_global_assessment, y = cox3_log)) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.05) +
  xlab("Physician Global Assessment") +
  ylab("COX3 (copies/uL)") +
  labs(
    title = "COX3 by Physician Global Assessment"
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_options
ggplot2::ggsave(paste0(output_dir, "cox3_vs_pga.png"),
  dpi = 300, width = 2000, height = 1500, units = "px"
)

# ====================================================================
# Total cfDNA, COX3 and ND2 Distribution Curves and Normality Checking
# ====================================================================
print("Creating distribution graphs...")

plot.distribution <- function(x, x_label, binwidth, density_scaling) {
  file_string <- paste0("distribution_", x, ".png")
  p <- df %>%
    ggplot(aes(x = .data[[x]])) +
    geom_histogram(binwidth = binwidth, color = "white", alpha = 0.5) +
    geom_density(
      size = 1, alpha = 0.5,
      aes(y = ..density.. * density_scaling)
    ) +
    xlab(x_label) +
    ylab("Count") +
    labs(title = paste0(x_label, " Distribution")) +
    theme_options
  print(p)
  ggplot2::ggsave(paste0(output_dir, file_string), p,
    dpi = 300, width = 2000, height = 1500, units = "px"
  )
}

plot.distribution(
  x = "total_cfdna", x_label = "Total cfDNA",
  binwidth = 0.05, density_scaling = 10
)
plot.distribution(
  x = "cox3_sapphire", x_label = "COX3",
  binwidth = 100, density_scaling = 10000
)
plot.distribution(
  x = "cox3_log", x_label = "Log COX3",
  binwidth = 0.2, density_scaling = 10
)
plot.distribution(
  x = "nd2_sapphire", x_label = "ND2",
  binwidth = 100, density_scaling = 10000
)
plot.distribution(
  x = "nd2_log", x_label = "Log ND2",
  binwidth = 0.2, density_scaling = 10
)

shapiro.test(df$total_cfdna)
shapiro.test(df$cox3_sapphire)
shapiro.test(df$cox3_log)
shapiro.test(df$nd2_sapphire)
shapiro.test(df$nd2_log)

print("Complete.")

# ====================================================================
# Correlation Analysis
# ====================================================================
# cor.test(df$plt_lab, df$cox3_log,
#    method=c("pearson", "kendall", "spearman"))
# Pearson's product-moment correlation
#
# data:  df$plt_lab and df$cox3_log
# t = 1.7857, df = 84, p-value = 0.07776
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.02150743  0.38741697
# sample estimates:
#       cor
# 0.1912394

# cor.test(df$haemoglobin, df$cox3_log,
#    method=c("pearson", "kendall", "spearman"))
# Pearson's product-moment correlation
#
# data:  df$haemoglobin and df$cox3_log
# t = -2.1724, df = 84, p-value = 0.03264
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.42189751 -0.01972763
# sample estimates:
#        cor
# -0.2306389

library("ggpubr")

plot.scatter <- function(df, x, y, x_label, y_label) {
  title_string <- paste0(y_label, " against ", x_label)
  file_string <- paste0(y, "_vs_", x)

  p <- ggscatter(df,
    x = x, y = y,
    add = "reg.line", conf.int = TRUE,
    cor.coef = TRUE, cor.method = "spearman", cor.coef.size = 8
  ) +
    xlab(x_label) +
    ylab(y_label) +
    labs(title = title_string) +
    theme_options
  print(p)
  ggplot2::ggsave(paste0(output_dir, file_string, ".png"), p,
    dpi = 300, width = 2500, height = 2000, units = "px"
  )
  return(p)
}

p_hb <- plot.scatter(df, x = "haemoglobin", y = "cox3_log",
  x_label = "Haemoglobin", y_label = "Log COX3")
p_plt <- plot.scatter(df, x = "plt_lab", y = "cox3_log",
  x_label = "Platelets", y_label = "Log COX3")
p_hb_plt <- gridExtra::grid.arrange(p_hb, p_plt, nrow = 1)
ggplot2::ggsave(paste0(output_dir, "hb_plt_correlates.png"), p_hb_plt,
  dpi = 300, width = 5000, height = 2000, units = "px"
)

plot.scatter(df, x = "plt_lab", y = "total_cfdna",
  x_label = "Platelets", y_label = "Total cfDNA")

plot.scatter(df, x = "haemoglobin", y = "total_cfdna",
  x_label = "Haemoglobin", y_label = "Total cfDNA")

p_calpro <- plot.scatter(df, x = "calprotectin", y = "cox3_log",
  x_label = "Calprotectin", y_label = "Log COX3")
p_crp <- plot.scatter(df, x = "crp", y = "cox3_log",
  x_label = "CRP", y_label = "Log COX3")
p_calpro_crp <- gridExtra::grid.arrange(p_calpro, p_crp, nrow = 1)
ggplot2::ggsave(paste0(output_dir, "known_biomarkers.png"), p_calpro_crp,
  dpi = 300, width = 5000, height = 2000, units = "px"
)

plot.scatter(df, x = "cucq_total", y = "cox3_log",
  x_label = "CUCQ-32", y_label = "Log COX3")

plot.scatter(df, x = "white_cell_count", y = "cox3_log",
  x_label = "White Cell Count", y_label = "Log COX3")

plot.scatter(df, x = "red_cell_count", y = "cox3_log",
  x_label = "Red Cell Count", y_label = "Log COX3")

plot.scatter(df, x = "age", y = "cox3_log",
  x_label = "Age", y_label = "Log COX3")

# ====================================================================
# Correlation Matrix
# ====================================================================
p1 <- df %>%
  select(c(
    "cox3_log",
    "total_cfdna",
    "nd2_log",
    "haemoglobin",
    "red_cell_count",
    "plt_lab",
    "calprotectin",
    "cucq_total"
  )) %>%
  ggpairs()
ggplot2::ggsave(paste0(output_dir, "correlation_matrix.png"), p1,
  dpi = 300, width = 4000, height = 4000, units = "px"
)

cor.test(df$nd2_sapphire, df$cox3_sapphire, method = "spearman")
cor.test(df$total_cfdna, df$cox3_sapphire, method = "spearman")

# ====================================================================
# Complete Mucosal Healing
# ====================================================================

t(prop.table(table(df_baseline$complete_mucosal_healing))) * 100
t(prop.table(table(df_baseline$endoscopic_improvement))) * 100


p_complete_mucosal_healing <- df_baseline %>%
  filter(!is.na(complete_mucosal_healing)) %>%
  ggplot(aes(x = complete_mucosal_healing, fill = complete_mucosal_healing)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100), width=0.5) +
  ylim(0,100) +
  xlab("Complete Mucosal Healing") +
  ylab("Percentage") +
  labs(
    title = "Complete Mucosal Healing at 3-6 months"
  ) +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  theme_options

p_endoscopic_improvement <- df_baseline %>%
  filter(!is.na(endoscopic_improvement)) %>%
  mutate(endoscopic_improvement = fct_relevel(endoscopic_improvement, 
                            "Yes", "No")) %>%
  ggplot(aes(x = endoscopic_improvement, fill = endoscopic_improvement)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100), width=0.5) +
  ylim(0,100) +
  xlab("Endoscopic Improvement") +
  ylab("Percentage") +
  labs(
    title = "Endoscopic Improvement at 3-6 months"
  ) +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  theme_options

p_complete_mucosal_healing
p_endoscopic_improvement

p_mucosal_healing <- gridExtra::grid.arrange(p_complete_mucosal_healing, p_endoscopic_improvement, nrow = 1)
ggplot2::ggsave(paste0(output_dir, "endoscopic_healing.png"), p_mucosal_healing,
                dpi = 300, width = 5000, height = 2000, units = "px"
)

df_baseline %>% 
  filter(!is.na(complete_mucosal_healing)) %>%
  group_by(study_group_name, complete_mucosal_healing) %>%
  count()

df_baseline %>% 
  filter(!is.na(complete_mucosal_healing)) %>%
  group_by(study_group_name, complete_mucosal_healing) %>%
  summarise(Percentage=n()) %>% 
  group_by(study_group_name) %>% 
  mutate(Percentage=Percentage/sum(Percentage)*100)

df_baseline %>% 
  filter(!is.na(endoscopic_improvement)) %>%
  group_by(study_group_name, endoscopic_improvement) %>%
  count()

df_baseline %>% 
  filter(!is.na(endoscopic_improvement)) %>%
  group_by(study_group_name, endoscopic_improvement) %>%
  summarise(Percentage=n()) %>% 
  group_by(study_group_name) %>% 
  mutate(Percentage=Percentage/sum(Percentage)*100)

print("Script successfully completed.")
