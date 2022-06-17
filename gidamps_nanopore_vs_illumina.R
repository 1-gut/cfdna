# GIDAMPs Nanopore vs Illumina
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
input_csv <- "data/gidamps/nanopore_vs_illumina/cleaned_nanopore_vs_illumina_data.csv"
gidamps_clinical_data <- "data/gidamps/clinical/merged_df.rds"
output_dir <- "output/gidamps/nanopore_vs_illumina/"
# -----------------------------------------------------------------------------
# Description
# -----------------------------------------------------------------------------
# Loads data from nanopore vs illumina sequencing performed in Apr 2021
# Plots raw reads, classified reads and microbial reads
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggsignif)
library(gridExtra)
source("./theme_options.R")

df <- read.csv(input_csv)
df_gidamps <- readRDS(gidamps_clinical_data)

df$sampling_date <- as.Date(df$sampling_date)

df <- left_join(df, df_gidamps, by = c("study_id", "sampling_date"))

df$raw_reads_per_million <- df$raw_reads / 1000000

# Plot raw reads
p_raw_reads_each_sample <- df %>%
  ggplot(aes(x = reorder(study_id, desc(raw_reads_per_million)),
    y = raw_reads_per_million, fill = group.x)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(label = comma) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    x = "Study ID",
    y = "Raw reads (million)",
    title = "Raw reads for each sample"
  ) +
  theme_options +
  theme(legend.position = "top")

t.test(raw_reads_per_million ~ group.x, data = df)

p_raw_reads_aggregated <- df %>%
  ggplot(aes(x = group.x, y = raw_reads_per_million, fill = group.x)) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(width = 0.3) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    x = "Sequencing Technology",
    y = "Raw reads (million)",
    title = "Raw reads aggregated"
  ) +
  theme_options +
  geom_signif(
    comparisons = list(c("Illumina", "Nanopore")),
    y_position = 7,
    annotations = "p < 0.001",
    tip_length = 0.01,
    vjust = 0
  )

p_raw_reads <- grid.arrange(
  p_raw_reads_each_sample,
  p_raw_reads_aggregated,
  nrow = 1
)
ggsave(paste0(output_dir, "raw_reads.png"), p_raw_reads,
  dpi = 300, width = 5000, height = 2000, units = "px")

# classification performance

p_classified_reads_each_sample <- df %>%
  ggplot(aes(x = reorder(study_id, desc(classified_reads)),
    y = classified_reads, fill = group.x)) +
  geom_bar(stat = "identity") +
  ylim(0, 100) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    x = "Study ID",
    y = "Classified reads (percentage)",
    title = "Classified reads for each sample"
  ) +
  theme_options +
  theme(legend.position = "top")
print(p_classified_reads_each_sample)

t.test(classified_reads ~ group.x, data = df)

p_classified_reads_aggregated <- df %>%
  ggplot(aes(x = group.x, y = classified_reads, fill = group.x)) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(width = 0.3) +
  ylim(0, 100) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    x = "Sequencing Technology",
    y = "Classified reads (percentage)",
    title = "Classified reads aggregated"
  ) +
  theme_options +
  geom_signif(
    comparisons = list(c("Illumina", "Nanopore")),
    y_position = 90,
    annotations = "p < 0.001",
    tip_length = 0.01,
    vjust = 0
  )
print(p_classified_reads_aggregated)
p_classified_reads <- grid.arrange(
  p_classified_reads_each_sample,
  p_classified_reads_aggregated,
  nrow = 1
)
ggsave(paste0(output_dir, "classified_reads.png"), p_classified_reads,
  dpi = 300, width = 5000, height = 2000, units = "px")

# Microbial reads
p_microbial_reads_each_sample <- df %>%
  ggplot(aes(x = reorder(study_id, desc(microbial_reads)),
    y = microbial_reads, fill = group.x)) +
  geom_bar(stat = "identity") +
  ylim(0, 3) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    x = "Study ID",
    y = "Microbial reads (percentage)",
    title = "Microbial reads for each sample"
  ) +
  theme_options +
  theme(legend.position = "top")
print(p_microbial_reads_each_sample)

t.test(microbial_reads ~ group.x, data = df)

p_microbial_reads_aggregated <- df %>%
  ggplot(aes(x = group.x, y = microbial_reads, fill = group.x)) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(width = 0.3) +
  ylim(0, 3) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    x = "Sequencing Technology",
    y = "Microbial reads (percentage)",
    title = "Microbial reads aggregated"
  ) +
  theme_options +
  geom_signif(
    comparisons = list(c("Illumina", "Nanopore")),
    y_position = 2.5,
    annotations = "p = 0.055",
    tip_length = 0.01,
    vjust = 0
  )
print(p_microbial_reads_aggregated)
p_microbial_reads <- grid.arrange(
  p_microbial_reads_each_sample,
  p_microbial_reads_aggregated,
  nrow = 1
)
ggsave(paste0(output_dir, "microbial_reads.png"), p_microbial_reads,
  dpi = 300, width = 5000, height = 2000, units = "px")



print("Script successfully completed.")
