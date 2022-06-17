output_dir <- "output/chr_counts/"
input_idxstats <- "data/03052022/reports_pipeline_v1.3/grch38/idxstats/" # against GRCh38
input_flagstat <- "data/03052022/reports_pipeline_v1.3/grch38/flagstat/" # against GRCh38
input_mito_flagstat <- "data/03052022/reports_pipeline_v1.3/mito/samtools_flagstat/" # against mito genome
sample_metadata <- read.csv("data/phyloseq_metadata.csv", fileEncoding = "UTF-8")
library(dplyr)
library(tibble)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(scales)
source("./theme_options.R")

# ------------------------------------------------------------------------
# Initial Dataframe Creation
# -------------------------------
# Goes through all the idxstats file and extracts the data along with
# using the filename as sample name
# then goes through all flagstat files to extract total mapped reads from
# line 7 and merges it all into a single df for downstream analysis
# ------------------------------------------------------------------------

tsv_files <- list.files(input_idxstats, full.names=TRUE)
flagstat_files <- list.files(input_flagstat, full.names=TRUE)
mito_flagstat_files <- list.files(input_mito_flagstat, full.names=TRUE)

# Make IDXSTATS Dataframe
chr_to_keep <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                 "chr22", "chrX", "chrY", "chrM")

df <- data.frame()

for (file in tsv_files) {
  file_name <- tail(stringr::str_split(file, "/")[[1]],1) # obtain full file name
  sample_name <- stringr::str_split(file_name, "_")[[1]][[1]] # get the sample name from the first part of the file name
  
  tsv_data <- read.table(file[[1]], sep="\t")
  colnames(tsv_data) <- c("reference_sequence_name", "sequence_length", "mapped_reads", "unmapped_reads")
  tsv_data <- tsv_data %>% filter(reference_sequence_name %in% chr_to_keep) # keep only desired chr
  tsv_data$mapped_reads_normalized <- tsv_data$mapped_reads / tsv_data$sequence_length * 1000000 # create normalized column to chr_length and scale by 1000000
  
  mapped_reads_df <- tsv_data %>%
    select(c("reference_sequence_name", "mapped_reads")) %>%
    t %>%
    as.data.frame
  colnames(mapped_reads_df) <- mapped_reads_df[1,]
  mapped_reads_df <- mapped_reads_df[2,]
  mapped_reads_df$sample_name <- sample_name
  
  mapped_reads_normalized_df <- tsv_data %>%
    select(c("reference_sequence_name", "mapped_reads_normalized")) %>%
    t %>%
    as.data.frame
  colnames(mapped_reads_normalized_df) <- mapped_reads_normalized_df[1,]
  mapped_reads_normalized_df <- mapped_reads_normalized_df[2,]
  mapped_reads_normalized_df$sample_name <- sample_name
  
  combined_df <- dplyr::full_join(mapped_reads_df, mapped_reads_normalized_df, by="sample_name", suffix=c("","_normalized"))
  
  df <- rbind(df, combined_df)
}

df <- df %>% relocate(sample_name, .before = chr1)

# add total mapped_reads from GRCh38

total_mapped_reads_df <- data.frame()

for (file in flagstat_files) {
  file_name <- tail(stringr::str_split(file, "/")[[1]],1)
  sample_name <- stringr::str_split(file_name, "_")[[1]][[1]]

  flagstat_data <- read.table(file, sep="\n")
  total_mapped_reads <- as.numeric(stringr::str_split(flagstat_data[7,], " ")[[1]][[1]])
  flagstat_df <- data.frame(sample_name = sample_name, total_mapped_reads = total_mapped_reads)
  
  total_mapped_reads_df <- rbind(total_mapped_reads_df, flagstat_df)
}

df <- dplyr::left_join(df, total_mapped_reads_df, by="sample_name")
df <- df %>% relocate(total_mapped_reads, .before = chr1)

# add metadata

sample_metadata <- sample_metadata %>%
  select(c("study_id", "study_group_name", "ibd_status_collapsed")) %>%
  dplyr::rename(sample_name = study_id)


df <- dplyr::left_join(df, sample_metadata, by="sample_name")
df$ibd_status_collapsed <- as.factor(df$ibd_status_collapsed)
df$ibd_status_collapsed <- relevel(df$ibd_status_collapsed, "Remission")
df$ibd_status_collapsed <- relevel(df$ibd_status_collapsed, "Active")

# add ChrM reads from mito pipeline

mito_df <- data.frame()
for (file in mito_flagstat_files) {
  file_name <- tail(stringr::str_split(file, "/")[[1]],1)
  sample_name <- stringr::str_split(file_name, "_")[[1]][[1]]
  
  flagstat_data <- read.table(file, sep="\n")
  chrM_mito <- as.numeric(stringr::str_split(flagstat_data[7,], " ")[[1]][[1]])
  mito_flagstat_df <- data.frame(sample_name = sample_name, chrM_mito = chrM_mito)
  
  mito_df <- rbind(mito_df, mito_flagstat_df)
}

df <- dplyr::left_join(df, mito_df, by="sample_name")

# ------------------------------------------------------------------------
# Begin Analysis Here
# ------------------------------------------------------------------------

# Go through every chromosome column and compute percentage
list_of_columns <- colnames(df)
list_of_columns <- list_of_columns[
  !list_of_columns %in% c("sample_name", "total_mapped_reads", "study_group_name", "ibd_status_collapsed")
]
for (column in list_of_columns) {
  df[[column]] <- as.numeric(df[[column]])
  df[[paste(column, "percent", sep = "_")]] <- paste0((df[[column]] / df$total_mapped_reads * 100))
  df[[paste(column, "percent", sep = "_")]] <- as.numeric(df[[paste(column, "percent", sep = "_")]])
}

# plot all host chromosomes

list_of_percentage_columns <- colnames(df)
list_of_percentage_columns <- list_of_percentage_columns[grepl("^.*\\_normalized_percent$", list_of_percentage_columns)]
# 
# percentage_plot_list <- list()
# 
# for (chr_column in list_of_percentage_columns) {
#   chr_column_label <- stringr::str_split(chr_column, "_")[[1]][[1]]
#   print(chr_column)
#   p <- df %>%
#     filter(sample_name != "NC") %>%
#     filter(sample_name != "PC") %>%
#     ggplot(aes_string(y = chr_column, x = "sample_name")) +
#     geom_bar(stat = "identity") +
#     ylim(0,4) +
#     labs(x="Sample ID", y=chr_column_label) +
#     theme_options
#   percentage_plot_list[[chr_column]] <- p
# }
# 
# g2 <- grid.arrange(grobs = percentage_plot_list, ncol = 4)

# ------------------------------------------------------------------------
# Visualize all chromosomes mapping ratio except chrM
# ------------------------------------------------------------------------

percentage_plot_list_without_chrm <- list()
list_of_percentage_columns_without_chrm <- head(list_of_percentage_columns, -1)
for (chr_column in list_of_percentage_columns_without_chrm) {
  chr_column_label <- stringr::str_split(chr_column, "_")[[1]][[1]]
  print(chr_column)
  p <- df %>%
    filter(sample_name != "NC") %>%
    filter(sample_name != "PC") %>%
    ggplot(aes_string(y = chr_column, x = "sample_name")) +
    geom_bar(stat = "identity") +
    ylim(0,0.05) +
    labs(x="Sample ID", y=chr_column_label)+
    theme_options
  percentage_plot_list_without_chrm[[chr_column]] <- p
}

g3 <- grid.arrange(grobs = percentage_plot_list_without_chrm, ncol = 4)
ggsave(paste0(output_dir, "all_chromosomes_normalized.png"), g3,
  width = 10000, height = 8000, units = "px", limitsize = FALSE
)

# ------------------------------------------------------------------------
# Visualize chrM ratio vs chr20
# ------------------------------------------------------------------------

p_chrm <- df %>%
  filter(sample_name != "NC") %>%
  filter(sample_name != "PC") %>%
  ggplot(aes_string(y = "chrM_normalized_percent", x = "sample_name")) +
  geom_bar(stat = "identity") +
  ylim(0,4) +
  labs(x="Sample ID", y="chrM (normalized ratio)") +
  theme_options

p_chr20 <- df %>%
  filter(sample_name != "NC") %>%
  filter(sample_name != "PC") %>%
  ggplot(aes_string(y = "chr20_normalized_percent", x = "sample_name")) +
  geom_bar(stat = "identity") +
  ylim(0,4) +
  labs(x="Sample ID", y="chr20 (normalized ratio)") +
  theme_options

g4 <- grid.arrange(p_chrm, p_chr20, nrow=1)
ggsave(paste0(output_dir, "chr20_vs_chrM.png"), g4, dpi = 300, width = 5000, height = 2000, units = "px")


# ------------------------------------------------------------------------
# Plot mitochondria reads
# ------------------------------------------------------------------------

plot_mito_activity <- df %>%
  filter(sample_name != c("NC")) %>%
  filter(sample_name != "PC") %>%
  group_by(ibd_status_collapsed) %>%
  ggplot(aes(y = chrM_percent, x = ibd_status_collapsed, fill=ibd_status_collapsed)) +
  geom_boxplot(width=0.5) +
  labs(x="IBD Activity", y="Mitochondrial Reads (%)", title="Mitochondrial Reads by IBD Activity") +
  scale_fill_brewer(palette="Pastel1") +
  ylim(0,0.2) +
  theme_options
# ggsave(paste0(output_dir, "mitochondrial_reads_by_activity.png"), plot_mito_activity, dpi = 300, width = 2500, height = 2000, units = "px")

# Plot mito reads per sample

plot_mito_percent <- df %>%
  subset(sample_name != c("NC")) %>%
  subset(sample_name != "PC") %>%
  ggplot(aes(x = reorder(sample_name, -total_mapped_reads), y = chrM_percent)) +
  geom_col() +
  labs(x="Sample ID", y="Mitochondrial Reads (%)", title="Percentage of mitochondrial reads per sample") +
  scale_fill_brewer(palette="Pastel1") +
  ylim(0,0.2) +
  theme_options
# ggsave(paste0(output_dir, "mitochondrial_reads_per_sample.png"), plot_mito_percent, dpi = 300, width = 2500, height = 2000, units = "px")


# ------------------------------------------------------------------------
# Plot bacterial reads
# ------------------------------------------------------------------------
df <- read.csv("data/03052022/reports_pipeline_v1.3/kraken2/k2report/pavian_visualization_raw_read_numbers.csv")


df$Bacterial.reads <- as.numeric(gsub(",", "", df$Bacterial.reads))
df$Number.of.raw.reads <- as.numeric(gsub(",", "", df$Number.of.raw.reads))
df$bacterial_percentage <- df$Bacterial.reads / df$Number.of.raw.reads * 100
df$Name <- sapply(df$Name, function(x){stringr::str_split(x, "_")[[1]][[1]]})

sample_metadata <- read.csv("data/phyloseq_metadata.csv", fileEncoding = "UTF-8")
sample_metadata <- sample_metadata %>%
  select(c("study_id", "study_group_name", "ibd_status_collapsed")) %>%
  dplyr::rename(Name = study_id)

df <- left_join(df, sample_metadata, by="Name")

df$ibd_status_collapsed <- as.factor(df$ibd_status_collapsed)
df$ibd_status_collapsed <- relevel(df$ibd_status_collapsed, "Remission")
df$ibd_status_collapsed <- relevel(df$ibd_status_collapsed, "Active")

plot_bacterial_percent <- df %>% 
  subset(Name != "PC") %>%
  subset(Name != "NC") %>%
  ggplot(aes(x=reorder(Name, -Number.of.raw.reads), y=bacterial_percentage)) +
  geom_col() +
  labs(x="Sample ID", y="Bacterial Reads (%)", title="Percentage of bacterial reads per sample") +
  ylim(0,0.2) +
  theme_options
# ggsave(paste0(output_dir, "raw_bacterial_percentage.png"), plot_bacterial_percent, dpi = 300, width = 2500, height = 2000, units = "px")

# Anova shows no significant differences
df_without_controls <- df %>% subset(Name != "PC") %>% subset(Name != "NC")
bacterial_aov <- aov(bacterial_percentage ~ ibd_status_collapsed, data=df_without_controls)
summary(bacterial_aov)
tukey <- TukeyHSD(bacterial_aov)
tukey


plot_bacterial_activity <- df %>% 
  subset(Name != "PC") %>%
  subset(Name != "NC") %>%
  group_by(ibd_status_collapsed) %>%
  ggplot(aes(x=ibd_status_collapsed, y=bacterial_percentage, fill=ibd_status_collapsed)) +
  geom_boxplot(width=0.5) +
  labs(x="IBD Activity", y="Bacterial Reads (%)", title="Bacterial Reads by IBD Activity") +
  scale_fill_brewer(palette="Pastel1") +
  ylim(0,0.2) +
  theme_options
# ggsave(paste0(output_dir, "bacterial_reads_by_activity.png"), plot_bacterial_activity, dpi = 300, width = 2500, height = 2000, units = "px")

# Same as above just different title
plot_bacterial_activity_rel <- df %>% 
  subset(Name != "PC") %>%
  subset(Name != "NC") %>%
  group_by(ibd_status_collapsed) %>%
  ggplot(aes(x=ibd_status_collapsed, y=bacterial_percentage, fill=ibd_status_collapsed)) +
  geom_boxplot(width=0.5) +
  labs(x="IBD Activity", y="Bacterial Reads (%)", title="Bacterial Reads (Relative)") +
  scale_fill_brewer(palette="Pastel1") +
  ylim(0,0.2) +
  theme_options

plot_bacterial_activity_raw <- df %>% 
  subset(Name != "PC") %>%
  subset(Name != "NC") %>%
  group_by(ibd_status_collapsed) %>%
  ggplot(aes(x=ibd_status_collapsed, y=Bacterial.reads, fill=ibd_status_collapsed)) +
  geom_boxplot(width=0.5) +
  labs(x="IBD Activity", y="Bacterial Reads (Raw count))", title="Bacterial Reads (Raw)") +
  scale_y_continuous(labels=comma) +
  scale_fill_brewer(palette="Pastel1") +
  theme_options

# ------------------------------------------------------------------------
# Create mito and bacterial composite graphs
# ------------------------------------------------------------------------

plot_bacterial_reads_composite <- grid.arrange(plot_bacterial_activity_raw, plot_bacterial_activity_rel, nrow=1)
ggsave(paste0(output_dir, "bacterial_reads_composite.png"), plot_bacterial_reads_composite, dpi = 300, width = 5000, height = 2000, units = "px")

plot_bacterial_mito_composite <- grid.arrange(plot_bacterial_percent, plot_mito_percent, plot_bacterial_activity, plot_mito_activity, nrow=2)
ggsave(paste0(output_dir, "bacterial_mito_composite.png"), plot_bacterial_mito_composite, dpi = 300, width = 5000, height = 4000, units = "px")

print("Script completed.")
