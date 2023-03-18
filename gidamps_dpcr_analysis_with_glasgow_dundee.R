# GIDAMPS DPCR ANALYSIS
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
output_dir <- "output/gidamps/dpcr_glasgowdundee/"
cohort_dir <- "output/gidamps/cohort_glasgowdundee/"

input_dpcr_file <- "data/gidamps/dpcr/dpcr_analysis_14032023.csv"
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
library("ggpubr")
source("./theme_options.R")
print("Complete.")

# Load dpcr_data
print("Loading sample_date from original cfDNA extraction excel...")
dpcr <- read.csv(input_dpcr_file, fileEncoding = "UTF-8-BOM")

# here's the magic
dpcr$study_id <- sub("GID-*", "", dpcr$sample_number)
dpcr$study_id <- sub("*-P", "", dpcr$study_id)
dpcr$study_id <- sub("*-HC", "", dpcr$study_id)
dpcr$sampling_date <- as.Date(dpcr$sample_date, format = "%d/%m/%Y")

# Add Log transformed columns
dpcr$total_cfdna_log <- log10(dpcr$total_cfdna)
dpcr$cox3_log <- log10(dpcr$cox3)
dpcr$nd2_log <- log10(dpcr$nd2)
print("Complete.")


# Read GID dataframe in and format it for merging
print("Read GI-DAMPs data in...")
gid <- readRDS(input_gidamps_merged_df)
gid$sampling_date <- as.Date(gid$sampling_date, format = "%Y-%m-%d")
print("Complete.")

# Merge GIDAMPs onto dPCR Data
print("Creating merged file...")
df <- left_join(dpcr, gid, by = c("study_id", "sampling_date"))
print("Complete.")

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
    "dpcr_date",
    "sampling_date",
    )) %>%
  select(-starts_with("drug_level_"))

demographic_df <- subset(demographic_df, group != "NC") # remove NC
demographic_df <- subset(demographic_df, group != "IBDU") # remove IBDU
demographic_df$group[demographic_df$group == "Non-IBD"] <- "HC" # collapse HC and Non-IBD

# Fill NA with 0
demographic_df[,c("aza","mp", "ada", "vedo", "uste", "tofa")][is.na(demographic_df[,c("aza","mp", "ada", "vedo", "uste", "tofa")])] <- 0

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
  count(previous_thiopurine)
cat("Previous Biologic by Group")
ibd_demo %>%
  group_by(group) %>%
  count(previous_biologic)

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

p.cox3nd2 <- grid.arrange(p.cox3, p.cox3_log, p.nd2, p.nd2_log, nrow=2)
ggsave(paste0(output_dir, "cox3_nd2_distrbutions_combined.png"), p.cox3nd2, dpi = 300, width = 4000, height = 3000, units = "px")
p3 <- grid.arrange(p1, p2, nrow=1)
ggsave(paste0(output_dir, "total_cfdna_combined.png"), p3, dpi = 300, width = 4000, height = 1500, units = "px")
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
  ggsave(paste0(output_dir, file_string), p, dpi = 300, width = 2000, height = 1500, units = "px")
  return(p)
}

plot.by.group(y="cox3", y_label="COX3 (copies/uL)", title="COX3 by Groups (Log scale)", log_scale=TRUE)
plot.by.group(y="nd2", y_label="ND2 (copies/uL)", title="ND2 by Groups (Log scale)", log_scale=TRUE)

p <- plot.by.group(y="total_cfdna", y_label="Total cfDNA (ng/uL)", title="Total cfDNA by Groups", log_scale=FALSE)
p + geom_signif(
  comparisons = list(c("Ulcerative colitis", "Healthy controls")),
  y_position = 3,
  annotations = "p < 0.05",
  tip_length = 0.01,
  vjust = 0
) 

# Remove negative controls significance bars
# +
#   geom_signif(
#     comparisons = list(c("Ulcerative colitis", "Negative controls")),
#     y_position = 2.3,
#     annotations = "***",
#     tip_length = 0,
#     vjust = 0.5
#   ) +
#   geom_signif(
#     comparisons = list(c("Crohn's disease", "Negative controls")),
#     y_position = 2.4,
#     annotations = "***",
#     tip_length = 0,
#     vjust = 0.5
#   ) +
#   geom_signif(
#     comparisons = list(c("Healthy controls", "Negative controls")),
#     y_position = 2.2,
#     annotations = "***",
#     tip_length = 0,
#     vjust = 0.5
#   )
ggsave(paste0(output_dir, "by_group_total_cfdna.png"), dpi = 300, width = 2500, height = 2000, units = "px")



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
# Remission     0.7748                -         -     
# Active        0.7136                0.5015    -     
# Highly active 0.0191                0.0177    0.0024

p <- plot.by.ibd.status(
  df=ibd_status_df_exclude_outlier_cfdna,
  y="total_cfdna",
  y_label="Total cfDNA (ng/uL)",
  title="Total cfDNA by Activity",
  caption="* p<0.05 ** p<0.01.\nKruskal-Wallis p<0.001.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied.",
  log_scale=FALSE
)
p +
  geom_signif(
    comparisons = list(c("Highly active", "Active")),
    y_position = 2.6,
    tip_length = 0,
    annotations = "**",
    vjust = 0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Remission")),
    y_position = 2.7,
    tip_length = 0,
    annotations = "*",
    vjust = 0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Biochemical remission")),
    annotations = "*",
    y_position = 2.8,
    tip_length = 0,
    vjust = 0.5
  )
ggsave(paste0(output_dir, "activity_cfdna.png"), dpi = 300, width = 2500, height = 2000, units = "px")


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

# Median and IQR for cox3 and total cfDNA

ibd_status_df_exclude_outlier_cfdna %>%
  group_by(ibd_status) %>%
  summarise_at(vars(cox3), list(median = median, iqr = IQR), na.rm=TRUE)

ibd_status_df_exclude_outlier_cfdna %>%
  group_by(ibd_status) %>%
  summarise_at(vars(nd2), list(median = median, iqr = IQR))

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


kruskal.test(nd2 ~ ibd_status, data = ibd_status_df_exclude_outlier_cfdna)
# Kruskal-Wallis rank sum test
# data:  nd2 by ibd_status
# Kruskal-Wallis chi-squared = 15.33, df = 3, p-value = 0.001555

pairwise.wilcox.test(
  ibd_status_df_exclude_outlier_cfdna$nd2,
  ibd_status_df_exclude_outlier_cfdna$ibd_status,
  p.adjust.method = "BH"
)
# Pairwise comparisons using Wilcoxon rank sum exact test 
# data:  ibd_status_df_exclude_outlier_cfdna$nd2 and ibd_status_df_exclude_outlier_cfdna$ibd_status 
#                 Biochemical remission Remission Active
#   Remission     0.731                 -         -     
#   Active        0.284                 0.347     -     
#   Highly active 0.044                 0.066     0.025 

p <- plot.by.ibd.status(
  df=ibd_status_df_exclude_outlier_cfdna,
  y="nd2",
  y_label="ND2 (copies/uL)",
  title="ND2 by Activity",
  caption="* p<0.05  ** p<0.01\nKruskal-Wallis p<0.05.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied.",
  log_scale=TRUE
)
p +
  geom_signif(
    comparisons = list(c("Highly active", "Active")),
    y_position = 4.1,
    tip_length = 0,
    annotations = "*",
    vjust=0.5
  ) +
  geom_signif(
    comparisons = list(c("Highly active", "Biochemical remission")),
    annotations = "*",
    y_position = 4.2,
    tip_length = 0,
    vjust=0.5
  )
ggsave(paste0(output_dir, "activity_nd2.png"), dpi = 300, width = 2500, height = 2000, units = "px")
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
    annotations = "*",
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
    caption = "* p<0.05. Kruskal-Wallis p<0.05.\nPost-hoc pairwise Wilcoxon test with Benjamini & Hochberg correction applied."
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
  xlab("IBD Activity") +
  ylab("Total cfDNA (ng/uL)") +
  labs(
    title = "Ulcerative Colitis",
    caption = "Kruskal-Wallis p<0.05.\nPost-hoc testing did not reach statistical significance due to small numbers."
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
  ggplot(aes(x = ibd_status, y = cox3)) +
  geom_boxplot() +
  scale_y_log10() +
  geom_jitter(width = 0.1) +
  facet_grid(group ~ .) +
  xlab("IBD Activity") +
  ylab("COX3 (copies/uL)") +
  ggtitle("COX3 by Group and Activity (Log Scale)")
# ggsave("output/dpcr/cox3_by_group_and_activity_log.png", dpi = 300, width = 2560, height = 1920, units = "px")


print("Complete.")
# ====================================================================
# ND2
# ====================================================================

print("Creating ND2 graphs...")

# ====================================================================
# ND2 by IBD group and activity (IBDU Excluded)
# ====================================================================

# Crohn's group comparison
kruskal.test(nd2 ~ ibd_status, data = cd_df)
pairwise.wilcox.test(cd_df$nd2, cd_df$ibd_status,
                     p.adjust.method = "BH"
)
#                 Biochemical remission Remission Active
#   Remission     0.495                 -         -     
#   Active        0.207                 1.000     -     
#   Highly active 0.002                 0.021     0.002 

# UC group comparison
kruskal.test(nd2 ~ ibd_status, data = uc_df)
pairwise.wilcox.test(uc_df$nd2, uc_df$ibd_status,
                     p.adjust.method = "BH"
)
#                 Biochemical remission Remission Active
#   Remission     1.00                  -         -     
#   Active        0.87                  0.87      -     
#   Highly active 0.87                  0.87      0.87  

ibd_status_df_exclude_outlier_cfdna_without_ibdu %>%
  ggplot(aes(x = ibd_status, y = nd2)) +
  geom_boxplot() +
  scale_y_log10() +
  geom_jitter(width = 0.1) +
  facet_grid(group ~ .) +
  xlab("IBD Activity") +
  ylab("ND2 (copies/uL)") +
  ggtitle("ND2 by Group and Activity (Log Scale)")
# ggsave("output/dpcr/nd2_by_group_and_activity_log.png", dpi = 300, width = 2560, height = 1920, units = "px")

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

plot.scatter <- function(df, x, y, x_label, y_label) {
  title_string <- paste0(y_label, " against ", x_label)
  file_string <- paste0(y, "_vs_", x)
  
  p <- ggscatter(df, x = x, y = y, 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "spearman", cor.coef.size=8) +
    xlab(x_label) +
    ylab(y_label) +
    labs(title = title_string) +
    theme_options
  print(p)
  ggsave(paste0(output_dir, file_string, ".png"), p, dpi = 300, width = 2500, height = 2000, units = "px")
  return(p)
}

p.hb <-plot.scatter(correlation_all, x="haemoglobin", y="cox3_log", x_label="Haemoglobin", y_label="Log COX3")
p.plt <- plot.scatter(correlation_all, x="plt_lab", y="cox3_log", x_label="Platelets", y_label="Log COX3")
p.hb_plt <- gridExtra::grid.arrange(p.hb, p.plt, nrow=1)
ggsave(paste0(output_dir, "hb_plt_correlates.png"), p.hb_plt, dpi = 300, width = 5000, height = 2000, units = "px")

p.calpro <- plot.scatter(correlation_all, x="calprotectin", y="cox3_log", x_label="Calprotectin", y_label="Log COX3")
p.crp <- plot.scatter(correlation_all, x="crp", y="cox3_log", x_label="CRP", y_label="Log COX3")
p.calpro_crp <- gridExtra::grid.arrange(p.calpro, p.crp, nrow=1)
ggsave(paste0(output_dir, "known_biomarkers.png"), p.calpro_crp, dpi = 300, width = 5000, height = 2000, units = "px")

print("Script successfully completed.")
