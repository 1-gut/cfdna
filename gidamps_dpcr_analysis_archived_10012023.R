# GIDAMPS DPCR ANALYSIS
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
output_dir <- "output/gidamps/dpcr/"
cohort_dir <- "output/gidamps/cohort/"

# input_dpcr_file <- "data/gidamps/dpcr/dpcr_analysis_17032022.csv"
# Reanalysis 10012023 here:
input_dpcr_file <- "data/gidamps/reanalysis_10012023/dpcr_analysis_09012023.csv"
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
library(stringr)
library(ggplot2)
library(ggsignif)
library(summarytools)
library(GGally)
library(RColorBrewer)
library(gridExtra)
library(compareGroups)
library(scales)
library("ggpubr")
source("./theme_options.R")
print("Complete.")

# Load dpcr_data
print("Loading sample_date from original cfDNA extraction excel...")
dpcr <- read.csv(input_dpcr_file, fileEncoding = "UTF-8-BOM")
dpcr$study_id <- as.numeric(str_extract(dpcr$sample_number, "[0-9]+"))
dpcr$sampling_date <- as.Date(dpcr$sample_date, format = "%d/%m/%Y")
dpcr$study_id[dpcr$study_id > 10000] <- NA

# Add Log transformed columns
dpcr$total_cfdna_log <- log10(dpcr$total_cfdna)
dpcr$cox3_log <- log10(dpcr$cox3)
dpcr$nd2_log <- log10(dpcr$nd2)
dpcr$gapdh_log <- log10(dpcr$gapdh)
print("Complete.")


# Read GID dataframe in and format it for merging
print("Read GI-DAMPs data in...")
gid <- readRDS(input_gidamps_merged_df)
gid$study_id <- as.numeric(gid$study_id) # BEWARE! This step drops all the Glasgow patients 136-etc.
gid$sampling_date <- as.Date(gid$sampling_date, format = "%Y-%m-%d")
print("Complete.")

# Merge GIDAMPs onto dPCR Data
print("Creating merged file...")
df <- left_join(dpcr, gid, by = c("study_id", "sampling_date"))
print("Complete.")

# Add correction columns
df$nd2_corrected_hb <- df$nd2 / df$haemoglobin
df$nd2_corrected_plt <- df$nd2 / df$plt_lab
df$nd2_corrected_gapdh <- df$nd2 / df$gapdh

# Label negative controls
levels(df$group) <- c(levels(df$group), "NC")
df$group[is.na(df$study_id)] <- "NC"

# ====================================================================
# Data Frame Creation (Remove awaiting diagnosis and exclude outliers)
# ====================================================================
df <- subset(df, group != "Awaiting diagnosis")
total_cfdna_exclude_outliers <- subset(df, total_cfdna < 10)

# Simplified DF (combine non-IBD and HC, remove IBDU from graphs)
simplified_df <- total_cfdna_exclude_outliers
simplified_df$group <- droplevels(simplified_df$group)
levels(simplified_df$group) <- c("CD", "UC", "IBDU", "HC", "HC", "NC")
simplified_df <- subset(simplified_df, group != "IBDU")
simplified_df$group <- droplevels(simplified_df$group)
levels(simplified_df$group) <- c("Crohn's disease", "Ulcerative colitis", "Healthy controls", "Negative controls")

# IBD DF
ibd_df <- subset(df, group == "UC" | group == "CD" | group == "IBDU")
# IBD DF where status NA is dropped
ibd_status_df <- subset(ibd_df, !is.na(ibd_status))
# Drop not applicable from dataframe
ibd_status_df$ibd_status <- droplevels(ibd_status_df$ibd_status)
# Exclude outlier (2 data points)
ibd_status_df_exclude_outlier_cfdna <- subset(ibd_status_df, total_cfdna < 10)
ibd_status_df_exclude_outlier_cfdna_without_ibdu <- subset(ibd_status_df_exclude_outlier_cfdna, group != "IBDU")

# ====================================================================
# Cohort Demographics
# ====================================================================
# First, let's clean up unused columns so output is easier to read
print("Saving demographic summaries...")
input_df <- total_cfdna_exclude_outliers
demographic_df <- input_df %>%
  select(-starts_with("tube_")) %>%
  select(-starts_with("redcap_")) %>%
  select(-c(
    "sampls_date_experiment4",
    "bloodtestdate_as_lbt",
    "biopsies_taken",
    "bloods_taken",
    "baseline_eims____1000",
    "baseline_mont_cd_beh____1000",
    "baseline_mont_cd_loc____1000",
    "baseline_recruitment_type",
    "data_entry_date",
    "questionnaire_comments",
    "hbinumber_of_liquid_stools",
    "hbiabdominal_mass1",
    "hbigeneral_well_being_as",
    "hbiabdominal_pain",
    "date_of_nhs_bloods",
    "adalimumab_test",
    "infliximab_test",
    "vedolizumab_test",
    "ustekinumab_test",
    "legacy_study_id",
    "diagnosis_note",
    "registration_location",
    "is_healthy_control",
    "medication_comments",
    "date_of_nhs_calprotectin",
    "baseline_gi_symptoms_desc",
    "stools_taken",
    "date_of_nhs_tdm",
    "faecal_test_date",
    "date_of_radiology",
    "ct_result",
    "mri_sb_result",
    "mri_pelvis_result",
    "endoscopy_result",
    "pathology_result",
    "radiology_comments",
    "consent_date",
    "past_medical_history_detail",
    "past_ibd_surgery_detail",
    "tonsillectomy_date",
    "appendicectomy_date",
    "sample_date",
    "dpcr_date"
    )) %>%
  select(-starts_with("drug_level_"))

demographic_df <- subset(demographic_df, group != "NC") # remove NC
demographic_df <- subset(demographic_df, group != "IBDU") # remove IBDU
demographic_df$group[demographic_df$group == "Non-IBD"] <- "HC" # collapse HC and Non-IBD

# Fill NA with 0
demographic_df[,c("aza","mp", "ada", "vedo", "uste", "tofa")][is.na(demographic_df[,c("aza","mp", "ada", "vedo", "uste", "tofa")])] <- 0
demographic_df[,c("iv_steroids","antibiotics", "rescue_therapy")][is.na(demographic_df[,c("iv_steroids","antibiotics", "rescue_therapy")])] <- 0

uc_demo <- subset(demographic_df, group == "UC")
cd_demo <- subset(demographic_df, group == "CD")
hc_demo <- subset(demographic_df, group == "HC")
ibd_demo <- subset(demographic_df, group == "UC" | group == "CD")

# ====================================================================
# Previous Thiopurine and Biologic Exposure Stats
# ====================================================================


# Create previous_thiopurine and previous_biologic columns
ibd_demo$previous_thiopurine <- ifelse(ibd_demo$aza == 1 | ibd_demo$mp == 1, 1, 0)
ibd_demo$previous_thiopurine[is.na(ibd_demo$previous_thiopurine)] <- 0
ibd_demo$previous_biologic <- ifelse(ibd_demo$ifx == "yes" | ibd_demo$ada == 1 | ibd_demo$vedo == 1 | ibd_demo$uste == 1 | ibd_demo$tofa == 1, 1, 0)
ibd_demo$previous_biologic[is.na(ibd_demo$previous_biologic)] <- 0

# Sense check
ibd_demo %>%
  select(c("aza", "mp", "previous_thiopurine"))
ibd_demo %>%
  select(c("ifx", "ada", "vedo", "uste", "tofa", "previous_biologic"))

# View previous_thiopurine
cat("Previous Thiopurine by Group")
ibd_demo %>%
  group_by(group) %>%
  dplyr::count(previous_thiopurine)
cat("Previous Biologic by Group")
ibd_demo %>%
  group_by(group) %>%
  dplyr::count(previous_biologic)

# ====================================================================
# Cohort comparison statistics
# ====================================================================

sink(file=paste0(cohort_dir, "cohort_statistics.txt"))
cat("GI-DAMPs Cohort Comparison Statistics\n")
cat("===================================================================\n")
cat("Previous Treatment Exposure\n")
cat("===================================================================\n")
previous.treatment.exposure <- compareGroups(group ~ previous_biologic + previous_thiopurine + ifx + vedo + uste + tofa + ada + aza + mp, data = ibd_demo)
createTable(previous.treatment.exposure)
cat("===================================================================\n")
cat("Clinical Parameters\n")
cat("===================================================================\n")
res <- compareGroups(group ~ haemoglobin + white_cell_count + plt_lab + albumin + crp + calprotectin, data = ibd_demo)
createTable(res)
cat("===================================================================\n")
cat("Drug Therapy at Sampling\n")
cat("===================================================================\n")
drug.therapy.at.sampling <- compareGroups(group ~ iv_steroids + rescue_therapy + antibiotics, data = ibd_demo)
createTable(drug.therapy.at.sampling)
cat("===================================================================\n")
cat("Overall IBD Activity Status\n")
cat("===================================================================\n")
ibd.status <- compareGroups(group ~ ibd_status, data = ibd_demo)
createTable(ibd.status)
cat("===================================================================\n")
cat("Age comparison between CD, UC and HC\n")
cat("===================================================================\n")
# 3-way comparisons
age.aov <- aov(age ~ group, data = demographic_df)
summary(age.aov)
pairwise.t.test(demographic_df$age, demographic_df$group,
                p.adjust.method = "BH")
cat("===================================================================\n")
cat("Sex comparison between CD, UC and HC\n")
cat("===================================================================\n")
chisq.test(demographic_df$sex,
           demographic_df$group)
sink()
# ====================================================================
# End cohort comparison statistics
# ====================================================================
# Summary Tools Visualization
view(dfSummary(demographic_df), file=paste0(cohort_dir, "dpcr_all.html"))
view(dfSummary(uc_demo), file=paste0(cohort_dir, "dpcr_uc.html"))
view(dfSummary(cd_demo), file=paste0(cohort_dir, "dpcr_cd.html"))
view(dfSummary(hc_demo), file=paste0(cohort_dir, "dpcr_hc.html"))

# Group Counts
demographic_df %>% count(group)

print("Completed.")
# ====================================================================
# Total cfDNA, COX3 and ND2 Distribution Curves and Normality Checking
# ====================================================================
print("Creating distribution graphs...")

plot.distribution <- function(df, x, x_label, binwidth, density_scaling, caption=NULL) {
  # shapiro <- shapiro.test(df[[x]])
  # capture.output(shapiro, file=paste0("output/gidamps/shapiro_", x, ".txt"))

  file_string <- paste0("distribution_", x, ".png")
  p <- df %>%
    ggplot(aes(x = .data[[x]])) +
    geom_histogram(binwidth = binwidth, color = "white", alpha = 0.5) +
    geom_density(size = 1, alpha = 0.5, aes(y = ..density.. * density_scaling)) +
    xlab(x_label) +
    ylab("Count") +
    labs(title = paste0(x_label, " Distribution"), caption=caption) +
    theme_options
  print(p)
  ggsave(paste0(output_dir, file_string), p, dpi = 300, width = 2000, height = 1500, units = "px")
  return(p)
}

p1 <- plot.distribution(df=total_cfdna_exclude_outliers, x="total_cfdna", x_label="Total cfDNA", binwidth=0.05, density_scaling = 10)
p2 <- plot.distribution(df=total_cfdna_exclude_outliers, x="total_cfdna_log", x_label="Total cfDNA (Log)", binwidth=0.05, density_scaling = 5)
p.cox3 <- plot.distribution(df=total_cfdna_exclude_outliers, x="cox3", x_label="COX3", binwidth=100, density_scaling = 10000)
p.cox3_log <- plot.distribution(df=total_cfdna_exclude_outliers, x="cox3_log", x_label="COX3 (Log)", binwidth=0.2, density_scaling = 10)
p.nd2 <- plot.distribution(df=total_cfdna_exclude_outliers, x="nd2", x_label="ND2", binwidth=100, density_scaling = 10000)
p.nd2_log <- plot.distribution(df=total_cfdna_exclude_outliers, x="nd2_log", x_label="ND2 (Log)", binwidth=0.2, density_scaling = 10)

p.gapdh <- plot.distribution(df=total_cfdna_exclude_outliers, x="gapdh", x_label="GAPDH", binwidth=10, density_scaling = 0.1)
p.gapdh_log <- plot.distribution(df=total_cfdna_exclude_outliers, x="gapdh_log", x_label="GAPDH (Log)", binwidth=0.1, density_scaling = 10)

p.cox3nd2 <- grid.arrange(p.cox3, p.cox3_log, p.nd2, p.nd2_log, nrow=2)
ggsave(paste0(output_dir, "cox3_nd2_distrbutions_combined.png"), p.cox3nd2, dpi = 300, width = 4000, height = 3000, units = "px")

p3 <- grid.arrange(p1, p2, nrow=1)
ggsave(paste0(output_dir, "total_cfdna_combined.png"), p3, dpi = 300, width = 4000, height = 1500, units = "px")

p.gapdh_combined <- grid.arrange(p.gapdh, p.gapdh_log, nrow=1)
ggsave(paste0(output_dir, "gapdh_combined.png"), p.gapdh_combined, dpi = 300, width = 4000, height = 1500, units = "px")
print("Complete.")

shapiro.test(total_cfdna_exclude_outliers$total_cfdna_log)

p4 <- total_cfdna_exclude_outliers %>% 
  ggplot(aes(sample=total_cfdna)) +
  stat_qq() +
  geom_qq() +
  geom_qq_line() +
  xlab("Theoretical") +
  ylab("Sample") +
  labs(title="Total cfDNA QQ Plot") +
  theme_options

p5 <- total_cfdna_exclude_outliers %>% 
  ggplot(aes(sample=total_cfdna_log)) +
  stat_qq() +
  geom_qq() +
  geom_qq_line() +
  xlab("Theoretical") +
  ylab("Sample") +
  labs(title="Total cfDNA (Log) QQ Plot") +
  theme_options

p6 <- grid.arrange(p1, p2, p4, p5, nrow=2)
ggsave(paste0(output_dir, "total_cfdna_combined_with_qq.png"), p6, dpi = 300, width = 4000, height = 3000, units = "px")

# ====================================================================
# Total cfDNA by Diagnosis Groups
# ====================================================================
print("Creating total cfDNA graphs...")
# ====================================================================
# Plot by Diagnosis Groups (Simplified DF)
# ====================================================================

plot.by.group <- function(y, y_label, title, caption=NULL, log_scale=FALSE) {
  file_string <- paste0("by_group_", y, ".png")
  p <-  simplified_df %>%
    ggplot(aes(x = group, y = .data[[y]], fill=group)) +
    geom_boxplot(outlier.shape=NA, width=0.5) +
    geom_jitter(width = 0.1, alpha=0.5) +
    xlab("Groups") +
    ylab(y_label) +
    labs(title=title, caption=caption) +
    scale_fill_brewer(palette="Pastel1") +
    theme_options 
  if (log_scale) {
    p <- p + scale_y_log10()
  }
  print(p)
  ggsave(paste0(output_dir, file_string), p, dpi = 300, width = 3000, height = 2000, units = "px")
  return(p)
}

plot.by.group(y="cox3", y_label="COX3 (copies/uL)", title="COX3 by Groups (Log scale)", log_scale=TRUE)
plot.by.group(y="nd2", y_label="ND2 (copies/uL)", title="ND2 by Groups (Log scale)", log_scale=TRUE)
plot.by.group(y="gapdh", y_label="GAPDH (copies/uL)", title="GAPDH by Groups (Log scale)", log_scale=TRUE)

p <- plot.by.group(
  y="total_cfdna",
  y_label="Total cfDNA (ng/uL)",
  title="Total cfDNA by Groups",
  caption="* p<0.05 ** p<0.01.\nKruskal-Wallis p<0.001.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied.",
  log_scale=FALSE)
p + geom_signif(
  comparisons = list(c("Ulcerative colitis", "Healthy controls")),
  y_position = 2.9,
  annotations = "**",
  tip_length = 0,
  vjust = 0.5
) + geom_signif(
  comparisons = list(c("Crohn's disease", "Healthy controls")),
  y_position = 3,
  annotations = "*",
  tip_length = 0,
  vjust = 0.5
) 
ggsave(paste0(output_dir, "by_group_total_cfdna.png"), dpi = 300, width = 3000, height = 2000, units = "px")



kruskal.test(total_cfdna ~ group, data = simplified_df)
pairwise.wilcox.test(
  simplified_df$total_cfdna,
  simplified_df$group,
  p.adjust.method = "BH"
)

kruskal.test(cox3 ~ group, data = simplified_df)
pairwise.wilcox.test(
  simplified_df$cox3,
  simplified_df$group,
  p.adjust.method = "BH"
)

kruskal.test(nd2 ~ group, data = simplified_df)
pairwise.wilcox.test(
  simplified_df$nd2,
  simplified_df$group,
  p.adjust.method = "BH"
)

kruskal.test(gapdh ~ group, data = simplified_df)
pairwise.wilcox.test(
  simplified_df$gapdh,
  simplified_df$group,
  p.adjust.method = "BH"
)


# ====================================================================
# Total cfDNA by IBD Status
# ====================================================================
plot.by.ibd.status <- function(df, y, y_label, title, caption=NULL, log_scale=FALSE) {
  
  p <-  df %>%
    ggplot(aes(x = ibd_status, y = .data[[y]], fill=ibd_status)) +
    geom_boxplot(outlier.shape=NA, width=0.5) +
    geom_jitter(width = 0.1, alpha=0.5) +
    xlab("IBD Activity") +
    ylab(y_label) +
    labs(title=title, caption=caption) +
    scale_fill_brewer(palette="Pastel1") +
    theme_options 
  
  if (log_scale) {
    p <- p + scale_y_log10()
  }
  print(p)
  return(p)
}



kruskal.test(total_cfdna ~ ibd_status, data = ibd_status_df_exclude_outlier_cfdna)
# Kruskal-Wallis rank sum test
# 
# data:  total_cfdna by ibd_status
# Kruskal-Wallis chi-squared = 21.314, df = 3, p-value = 9.059e-05
pairwise.wilcox.test(ibd_status_df_exclude_outlier_cfdna$total_cfdna, ibd_status_df_exclude_outlier_cfdna$ibd_status,
                     p.adjust.method = "BH"
)
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  ibd_status_df_exclude_outlier_cfdna$total_cfdna and ibd_status_df_exclude_outlier_cfdna$ibd_status 
# 
# Biochemical remission Remission Active 
# Remission     0.44059               -         -      
#   Active        0.79369               0.31398   -      
#   Highly active 0.00221               0.00089   3.5e-06

p <- plot.by.ibd.status(
  df=ibd_status_df_exclude_outlier_cfdna,
  y="total_cfdna",
  y_label="Total cfDNA (ng/uL)",
  title="Total cfDNA by Activity",
  caption="** p<0.01 *** p<0.001.\nKruskal-Wallis p<0.001.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied.",
  log_scale=FALSE
)
p +
  geom_signif(
    comparisons = list(c("Highly active", "Active")),
    y_position = 2.6,
    tip_length = 0,
    annotations = "***",
    vjust = 0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Remission")),
    y_position = 2.7,
    tip_length = 0,
    annotations = "***",
    vjust = 0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Biochemical remission")),
    annotations = "**",
    y_position = 2.8,
    tip_length = 0,
    vjust = 0.5
  )
ggsave(paste0(output_dir, "activity_cfdna.png"), dpi = 300, width = 3000, height = 2000, units = "px")

# ====================================================================
# COX3 by activity
# ====================================================================

kruskal.test(cox3 ~ ibd_status, data = ibd_status_df_exclude_outlier_cfdna)
# Kruskal-Wallis rank sum test
# 
# data:  cox3 by ibd_status
# Kruskal-Wallis chi-squared = 18.518, df = 3, p-value = 0.0003439

pairwise.wilcox.test(
  ibd_status_df_exclude_outlier_cfdna$cox3,
  ibd_status_df_exclude_outlier_cfdna$ibd_status,
  p.adjust.method = "BH"
)
# Pairwise comparisons using Wilcoxon rank sum exact test 
# 
# data:  ibd_status_df_exclude_outlier_cfdna$cox3 and ibd_status_df_exclude_outlier_cfdna$ibd_status 
# 
# Biochemical remission Remission Active
# Remission     0.5338                -         -     
#   Active        0.4173                0.4173    -     
#   Highly active 0.0343                0.0343    0.0051
# 
# P value adjustment method: B

# Median and IQR for nd2 and total cfDNA

ibd_status_df_exclude_outlier_cfdna %>%
  group_by(ibd_status) %>%
  summarise_at(vars(nd2), list(median = median, iqr = IQR))

ibd_status_df_exclude_outlier_cfdna %>%
  group_by(ibd_status) %>%
  summarise_at(vars(cox3), list(median = median, iqr = IQR), na.rm=TRUE)

ibd_status_df_exclude_outlier_cfdna %>%
  group_by(ibd_status) %>%
  summarise_at(vars(total_cfdna), list(median = median, iqr = IQR))



p <- plot.by.ibd.status(
  df=ibd_status_df_exclude_outlier_cfdna,
  y="cox3",
  y_label="COX3 (copies/uL)",
  title="COX3 by Activity",
  caption="* p<0.05 ** p<0.01.\nKruskal-Wallis p<0.001.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied.",
  log_scale=TRUE
)
p + geom_signif(
  comparisons = list(c("Highly active", "Active")),
  y_position = 4.2,
  tip_length = 0,
  annotations = "**",
  vjust = 0.5
) +
  geom_signif(
    comparisons = list(c("Highly active", "Remission")),
    y_position = 4.3,
    tip_length = 0,
    annotations = "*",
    vjust = 0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Biochemical remission")),
    annotations = "*",
    y_position = 4.4,
    tip_length = 0,
    vjust = 0.5
  )
ggsave(paste0(output_dir, "activity_cox3.png"), dpi = 300, width = 2500, height = 2000, units = "px")

# ====================================================================
# ND2 by activity
# ====================================================================

kruskal.test(nd2 ~ ibd_status, data = ibd_status_df_exclude_outlier_cfdna)

pairwise.wilcox.test(
  ibd_status_df_exclude_outlier_cfdna$nd2,
  ibd_status_df_exclude_outlier_cfdna$ibd_status,
  p.adjust.method = "BH"
)
# Pairwise comparisons using Wilcoxon rank sum exact test 
# 
# data:  ibd_status_df_exclude_outlier_cfdna$nd2 and ibd_status_df_exclude_outlier_cfdna$ibd_status 
# 
#                 Biochemical remission Remission Active
#   Remission     0.9215                -         -     
#   Active        0.3913                0.2206    -     
#   Highly active 0.0019                0.0018    0.0012

p <- plot.by.ibd.status(
  df=ibd_status_df_exclude_outlier_cfdna,
  y="nd2",
  y_label="ND2 (copies/uL)",
  title="ND2 by Activity",
  caption="** p<0.01\nKruskal-Wallis p<0.05.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied.",
  log_scale=TRUE
)
p +
  geom_signif(
    comparisons = list(c("Highly active", "Active")),
    y_position = 4.1,
    tip_length = 0,
    annotations = "**",
    vjust=0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Biochemical remission")),
    annotations = "**",
    y_position = 4.3,
    tip_length = 0,
    vjust=0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Remission")),
    annotations = "**",
    y_position = 4.2,
    tip_length = 0,
    vjust=0.5
  )
ggsave(paste0(output_dir, "activity_nd2.png"), dpi = 300, width = 3000, height = 2000, units = "px")
# ====================================================================
# GAPDH by activity
# ====================================================================

kruskal.test(gapdh_log ~ ibd_status, data = ibd_status_df_exclude_outlier_cfdna)

pairwise.wilcox.test(
  ibd_status_df_exclude_outlier_cfdna$gapdh,
  ibd_status_df_exclude_outlier_cfdna$ibd_status,
  p.adjust.method = "BH"
)

ibd_status_df_exclude_outlier_cfdna %>%
  group_by(ibd_status) %>%
  summarise_at(vars(gapdh), list(median = median, iqr = IQR), na.rm=TRUE)

p <- plot.by.ibd.status(
  df=ibd_status_df_exclude_outlier_cfdna,
  y="gapdh",
  y_label="GAPDH (copies/uL)",
  title="GAPDH by Activity",
  caption="* p<0.05\nKruskal-Wallis p<0.01.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied.",
  log_scale=TRUE
)
p +
  geom_signif(
    comparisons = list(c("Highly active", "Active")),
    y_position = 2.55,
    tip_length = 0,
    annotations = "*",
    vjust=0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Remission")),
    annotations = "*",
    y_position = 2.7,
    tip_length = 0,
    vjust=0.5
  )
ggsave(paste0(output_dir, "activity_gapdh.png"), dpi = 300, width = 3000, height = 2000, units = "px")
# ====================================================================
# ND2 Corrected Hb by activity
# ====================================================================
options(scipen=5)
kruskal.test(nd2_corrected_hb ~ ibd_status, data = ibd_status_df_exclude_outlier_cfdna)

pairwise.wilcox.test(
  ibd_status_df_exclude_outlier_cfdna$nd2_corrected_hb,
  ibd_status_df_exclude_outlier_cfdna$ibd_status,
  p.adjust.method = "BH"
)

p.nd2_hb <- plot.by.ibd.status(
  df=ibd_status_df_exclude_outlier_cfdna,
  y="nd2_corrected_hb",
  y_label="ND2 (copies/uL) divided by Hb (g/L)",
  title="ND2 by Activity Corrected for Haemoglobin",
  caption="* p<0.05 *** p<0.001\nKruskal-Wallis p<0.001.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied.",
  log_scale=TRUE
)
p.nd2_hb <- p.nd2_hb +
  geom_signif(
    comparisons = list(c("Highly active", "Active")),
    y_position = 2.1,
    tip_length = 0,
    annotations = "***",
    vjust=0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Biochemical remission")),
    annotations = "*",
    y_position = 2.4,
    tip_length = 0,
    vjust=0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Remission")),
    annotations = "***",
    y_position = 2.25,
    tip_length = 0,
    vjust=0.5
  ) +
  scale_y_log10(labels = label_comma())
ggsave(paste0(output_dir, "activity_nd2_corrected_hb.png"), dpi = 300, width = 3000, height = 2000, units = "px")

# ====================================================================
# ND2 Corrected plt by activity
# ====================================================================

kruskal.test(nd2_corrected_plt ~ ibd_status, data = ibd_status_df_exclude_outlier_cfdna)

pairwise.wilcox.test(
  ibd_status_df_exclude_outlier_cfdna$nd2_corrected_plt,
  ibd_status_df_exclude_outlier_cfdna$ibd_status,
  p.adjust.method = "BH"
)

p.nd2_plt <- plot.by.ibd.status(
  df=ibd_status_df_exclude_outlier_cfdna,
  y="nd2_corrected_plt",
  y_label="ND2 (copies/uL) divided by platelets",
  title="ND2 by Activity Corrected for Platelets",
  caption="* p<0.05\nKruskal-Wallis p<0.01.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied.",
  log_scale=TRUE
)
p.nd2_plt <- p.nd2_plt +
  geom_signif(
    comparisons = list(c("Highly active", "Active")),
    y_position = 2.1,
    tip_length = 0,
    annotations = "*",
    vjust=0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Remission")),
    annotations = "*",
    y_position = 2.25,
    tip_length = 0,
    vjust=0.5
  )
ggsave(paste0(output_dir, "activity_nd2_corrected_plt.png"), dpi = 300, width = 3000, height = 2000, units = "px")

p3 <- grid.arrange(p.nd2_hb, p.nd2_plt, nrow=2)
ggsave(paste0(output_dir, "activity_nd2_corrected_hb_plt.png"), p3, dpi = 300, width = 3000, height = 4000, units = "px")


# ====================================================================
# ND2 Corrected gapdh by activity
# ====================================================================

kruskal.test(nd2_corrected_gapdh ~ ibd_status, data = ibd_status_df_exclude_outlier_cfdna)

pairwise.wilcox.test(
  ibd_status_df_exclude_outlier_cfdna$nd2_corrected_gapdh,
  ibd_status_df_exclude_outlier_cfdna$ibd_status,
  p.adjust.method = "BH"
)

p <- plot.by.ibd.status(
  df=ibd_status_df_exclude_outlier_cfdna,
  y="nd2_corrected_gapdh",
  y_label="ND2 divided by GAPDH",
  title="ND2 by Activity Corrected for GAPDH",
  caption="Kruskal-Wallis no significant differences.",
  log_scale=TRUE
)
ggsave(paste0(output_dir, "activity_nd2_corrected_gapdh.png"), dpi = 300, width = 3000, height = 2000, units = "px")
# ====================================================================
# Total cfDNA by Group and Activity
# ====================================================================
# split df into UC and CD for statistics
cd_df <- subset(ibd_status_df_exclude_outlier_cfdna_without_ibdu, group == "CD")
uc_df <- subset(ibd_status_df_exclude_outlier_cfdna_without_ibdu, group == "UC")

p <- plot.by.ibd.status(
  df=cd_df,
  y="nd2",
  y_label="ND2 (copies/uL)",
  title="Crohn's Disease",
  log_scale=TRUE
)
p <- plot.by.ibd.status(
  df=uc_df,
  y="nd2",
  y_label="ND2 (copies/uL)",
  title="Ulcerative Colitis",
  log_scale=TRUE
)
p <- plot.by.ibd.status(
  df=cd_df,
  y="cox3",
  y_label="COX3 (copies/uL)",
  title="Crohn's Disease",
  log_scale=TRUE
)
p <- plot.by.ibd.status(
  df=uc_df,
  y="cox3",
  y_label="COX3 (copies/uL)",
  title="Ulcerative Colitis",
  log_scale=TRUE
)

# Crohn's group comparison
kruskal.test(total_cfdna ~ ibd_status, data = cd_df)
# Kruskal-Wallis rank sum test
# 
# data:  total_cfdna by ibd_status
# Kruskal-Wallis chi-squared = 11.799, df = 3, p-value = 0.008106

pairwise.wilcox.test(cd_df$total_cfdna, cd_df$ibd_status,
  p.adjust.method = "BH"
)
# Pairwise comparisons using Wilcoxon rank sum exact test 
# 
# data:  cd_df$total_cfdna and cd_df$ibd_status 
# 
# Biochemical remission Remission Active
# Remission     0.775                 -         -     
#   Active        0.775                 0.775     -     
#   Highly active 0.158                 0.031     0.024 
# 
# P value adjustment method: BH 

p1 <- ibd_status_df_exclude_outlier_cfdna %>%
  subset(group=="CD") %>%
  ggplot(aes(x = ibd_status, y = total_cfdna, fill=ibd_status)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width = 0.1, alpha=0.5) +
  geom_signif(
    comparisons = list(c("Highly active", "Active")),
    y_position = 2.6,
    tip_length = 0,
    annotations = "***",
    vjust = 0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Remission")),
    y_position = 2.7,
    tip_length = 0,
    annotations = "*",
    vjust = 0.5
  ) +
  xlab("IBD Activity") +
  ylab("Total cfDNA (ng/uL)") +
  labs(
    title = "Crohn's Disease",
    caption = "* p<0.05. *** p<0.001.\nKruskal-Wallis p<0.001.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied."
  ) +
  scale_fill_brewer(palette="Pastel1") +
  theme_options


# UC group comparison
kruskal.test(total_cfdna ~ ibd_status, data = uc_df)
# Kruskal-Wallis rank sum test
# data:  total_cfdna by ibd_status
# Kruskal-Wallis chi-squared = 8.3223, df = 3, p-value = 0.0398
pairwise.wilcox.test(uc_df$total_cfdna, uc_df$ibd_status,p.adjust.method = "BH")
# Pairwise comparisons using Wilcoxon rank sum exact test 
# 
# data:  uc_df$total_cfdna and uc_df$ibd_status 
# 
# Biochemical remission Remission Active
# Remission     1.00                  -         -     
#   Active        0.40                  0.65      -     
#   Highly active 0.17                  0.40      0.17  
# 
# P value adjustment method: BH 

p2 <- ibd_status_df_exclude_outlier_cfdna %>%
  subset(group=="UC") %>%
  ggplot(aes(x = ibd_status, y = total_cfdna, fill=ibd_status)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width = 0.1, alpha=0.5) +
  geom_signif(
    comparisons = list(c("Highly active", "Active")),
    y_position = 2.6,
    tip_length = 0,
    annotations = "*",
    vjust = 0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Biochemical remission")),
    y_position = 2.7,
    tip_length = 0,
    annotations = "*",
    vjust = 0.5
  ) +
  xlab("IBD Activity") +
  ylab("Total cfDNA (ng/uL)") +
  labs(
    title = "Ulcerative Colitis",
    caption = "* p<0.05.\nKruskal-Wallis p<0.001.\nPost-hoc testing did not reach statistical significance due to small numbers."
  ) +
  scale_fill_brewer(palette="Pastel1") +
  theme_options

p3 <- grid.arrange(p1, p2, nrow=1)
ggsave(paste0(output_dir, "total_cfdna_by_group_and_activity.png"), p3, dpi = 300, width = 5000, height = 2000, units = "px")

print("Complete.")
# ====================================================================
# COX3/ND2
# ====================================================================
print("Creating COX3/ND2 graphs...")

# ====================================================================
# ND2 by Group and Activity
# ====================================================================

# Crohn's group comparison
kruskal.test(nd2 ~ ibd_status, data = cd_df)
pairwise.wilcox.test(cd_df$nd2, cd_df$ibd_status,
                     p.adjust.method = "BH"
)
#                 Biochemical remission Remission Active 
#   Remission     0.92929               -         -      
#   Active        0.92929               0.92929   -      
#   Highly active 0.00033               0.00033   1.3e-05


p1 <- ibd_status_df_exclude_outlier_cfdna %>%
  subset(group=="CD") %>%
  ggplot(aes(x = ibd_status, y = nd2, fill=ibd_status)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width = 0.1, alpha=0.5) +
  geom_signif(
    comparisons = list(c("Highly active", "Active")),
    y_position = 4.2,
    tip_length = 0,
    annotations = "***",
    vjust = 0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Remission")),
    y_position = 4.35,
    tip_length = 0,
    annotations = "***",
    vjust = 0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Biochemical remission")),
    y_position = 4.5,
    tip_length = 0,
    annotations = "***",
    vjust = 0.5
  ) +
  xlab("IBD Activity") +
  ylab("ND2 (copies/uL)") +
  labs(
    title = "ND2 (Crohn's Disease)",
    caption = "*** p<0.001.\nKruskal-Wallis p<0.001.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied."
  ) +
  scale_y_log10() +
  scale_fill_brewer(palette="Pastel1") +
  theme_options

# UC group comparison
kruskal.test(nd2 ~ ibd_status, data = uc_df)
# p = 0.08873
pairwise.wilcox.test(uc_df$nd2, uc_df$ibd_status,p.adjust.method = "BH")
#                 Biochemical remission Remission Active
#   Remission     0.31                  -         -     
#   Active        0.17                  0.31      -     
#   Highly active 0.17                  0.31      0.31  

p2 <- ibd_status_df_exclude_outlier_cfdna %>%
  subset(group=="UC") %>%
  ggplot(aes(x = ibd_status, y = nd2, fill=ibd_status)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width = 0.1, alpha=0.5) +
  xlab("IBD Activity") +
  ylab("ND2 (copies/uL)") +
  labs(
    title = "ND2 (Ulcerative Colitis)",
    caption = "Kruskal-Wallis p=0.08.\nInter-group difference did not reach statistical significance."
  ) +
  scale_y_log10() +
  scale_fill_brewer(palette="Pastel1") +
  theme_options
p2
p3 <- grid.arrange(p1, p2, nrow=2)
ggsave(paste0(output_dir, "nd2_by_group_and_activity.png"), p3, dpi = 300, width = 3000, height = 4000, units = "px")


# ====================================================================
# COX3 by IBD group and activity (IBDU Excluded)
# ====================================================================

# Crohn's group comparison
kruskal.test(cox3 ~ ibd_status, data = cd_df)
pairwise.wilcox.test(cd_df$cox3, cd_df$ibd_status,
                     p.adjust.method = "BH"
)
#                 Biochemical remission Remission Active
#   Remission     0.229                 -         -     
#   Active        0.207                 1.000     -     
#   Highly active 0.002                 0.014     0.002 

# UC group comparison
kruskal.test(cox3 ~ ibd_status, data = uc_df)
pairwise.wilcox.test(uc_df$cox3, uc_df$ibd_status,
                     p.adjust.method = "BH"
)
#                 Biochemical remission Remission Active
#   Remission     1.00                  -         -     
#   Active        0.97                  0.51      -     
#   Highly active 0.93                  0.51      0.51 


# Log scale
ibd_status_df_exclude_outlier_cfdna_without_ibdu %>%
  ggplot(aes(x = ibd_status, y = cox3, fill=ibd_status)) +
  geom_boxplot() +
  scale_y_log10() +
  geom_jitter(width = 0.1, alpha=0.5) +
  facet_grid(. ~ group) +
  xlab("IBD Activity") +
  ylab("COX3 (copies/uL)") +
  ggtitle("COX3 by Group and Activity (Log Scale)") +
  scale_fill_brewer(palette="Pastel1") +
  theme_options
# ggsave("output/dpcr/cox3_by_group_and_activity_log.png", dpi = 300, width = 2560, height = 1920, units = "px")

print("Complete.")
# ====================================================================
# Correlation Analysis
# ====================================================================

# Start by looking at everybody
correlation_all <- subset(total_cfdna_exclude_outliers, group != "NC")

correlation_all %>%
  ggplot(aes(x = cox3, y = total_cfdna)) +
  geom_point()

# CD only Correlation
correlation_all$cox3_log <- log(correlation_all$cox3)
correlation_all$nd2_log <- log(correlation_all$nd2)

# Have explored correlation of cox3, nd2 and total_cfdna against most blood parameters
# No correlation found.
cor.test(correlation_all$total_cfdna, correlation_all$cox3, method = "spearman")
cor.test(correlation_all$nd2, correlation_all$cox3, method = "spearman")
cor.test(correlation_all$nd2, correlation_all$total_cfdna, method = "spearman")
cor.test(correlation_all$nd2, correlation_all$gapdh, method = "spearman")
cor.test(correlation_all$gapdh, correlation_all$total_cfdna, method = "spearman")


plot.scatter <- function(df, x, y, x_label, y_label, title_string="") {
  if (title_string=="") {
    title_string <- paste0(y_label, " against ", x_label)
  } else {
    title_string <- title_string
  }
  
  file_string <- paste0(y, "_vs_", x)
  
  p <- ggscatter(df, x = x, y = y, 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = FALSE, cor.method = "spearman", cor.coef.size=8) +
    stat_cor(
      method="spearman",
      p.accuracy = 0.001,
      r.digits = 2,
      size = 8
    ) +
    xlab(x_label) +
    ylab(y_label) +
    labs(title = title_string) +
    theme_options
  print(p)
  ggsave(paste0(output_dir, file_string, ".png"), p, dpi = 300, width = 2000, height = 1500, units = "px")
  return(p)
}

p.hb <-plot.scatter(correlation_all, x="haemoglobin", y="nd2_log", x_label="Haemoglobin", y_label="Log ND2")
p.plt <- plot.scatter(correlation_all, x="plt_lab", y="nd2_log", x_label="Platelets", y_label="Log ND2")
p.hb_plt <- gridExtra::grid.arrange(p.hb, p.plt, nrow=1)
ggsave(paste0(output_dir, "hb_plt_correlates.png"), p.hb_plt, dpi = 300, width = 4000, height = 1500, units = "px")

p.calpro <- plot.scatter(correlation_all, x="calprotectin", y="nd2_log", x_label="Calprotectin", y_label="Log ND2")
p.crp <- plot.scatter(correlation_all, x="crp", y="nd2_log", x_label="CRP", y_label="Log ND2")
p.crp <- ggpar(p.crp, xlim=c(0,100))
p.calpro_crp <- gridExtra::grid.arrange(p.calpro, p.crp, nrow=1)
ggsave(paste0(output_dir, "known_biomarkers.png"), p.calpro_crp, dpi = 300, width = 4000, height = 1500, units = "px")

p.nd2_gapdh <- plot.scatter(correlation_all, x="nd2_log", y="gapdh_log", x_label="Log ND2", y_label="Log GAPDH", title_string="GAPDH against ND2")
p.total_cfdna_gapdh <- plot.scatter(correlation_all, x="total_cfdna_log", y="gapdh_log", x_label="Log Total cfDNA", y_label="Log GAPDH", title_string="GAPDH against Total cfDNA")
p.gapdh_correlates <- gridExtra::grid.arrange(p.total_cfdna_gapdh, p.nd2_gapdh, nrow=1)
ggsave(paste0(output_dir, "gapdh_correlates.png"), p.gapdh_correlates, dpi = 300, width = 4000, height = 1500, units = "px")
print("Script successfully completed.")

