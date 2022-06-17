# COMBINED MUSIC GIDAMPS ANALYSIS FOR DPCR MEASUREMENT STATS
# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
output_dir <- "output/combined/"
# -----------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggprism)
library(ggplot2)
library(gridExtra)
library(ggpubr)

gidamps_dpcr <- read.csv("data/gidamps/dpcr/dpcr_analysis_17032022.csv", fileEncoding = "UTF-8-BOM")
music_dpcr <- read_excel("data/music/dpcr/music_cfdna_master_02062022.xlsx", sheet="cleaned_data")

theme_options <- theme_prism(base_size = 16) + theme(
  axis.title=element_text(size=21),
  plot.title=element_text(size=23),
  legend.position="none"
)

music_dpcr <- music_dpcr %>% rename(cox3 = cox3_sapphire, nd2 = nd2_sapphire)

df1 <- gidamps_dpcr %>% select(c("cox3","nd2", "total_cfdna"))
df2 <- music_dpcr %>% select(c("cox3","nd2", "total_cfdna"))

merged_df <- bind_rows(df1,df2)
merged_df <- subset(merged_df, total_cfdna < 10) # Exclude outliers

# Add log columns
merged_df$cox3_log <- log10(merged_df$cox3)
merged_df$nd2_log <- log10(merged_df$nd2)
merged_df$total_cfdna_log <- log10(merged_df$total_cfdna)

plot.distribution <- function(x, x_label, binwidth, density_scaling) {
  file_string <- paste0("distribution_", x, ".png")
  p <- merged_df %>%
    ggplot(aes(x = .data[[x]])) +
    geom_histogram(binwidth = binwidth, color = "white", alpha = 0.5) +
    geom_density(size = 1, alpha = 0.5, aes(y = ..density.. * density_scaling)) +
    xlab(x_label) +
    ylab("Count") +
    labs(title = paste0(x_label, " Distribution")) +
    theme_options
  print(p)
  ggsave(paste0(output_dir, file_string), p, dpi = 300, width = 2000, height = 1500, units = "px")
  return(p)
}

p.cox3 <- plot.distribution(x="cox3", x_label="COX3", binwidth=100, density_scaling=10000)
p.nd2 <- plot.distribution(x="nd2", x_label="ND2", binwidth=100, density_scaling=10000)
p.cox3_log <- plot.distribution(x="cox3_log", x_label="COX3 (Log)", binwidth=0.2, density_scaling=10)
p.nd2_log <- plot.distribution(x="nd2_log", x_label="ND2 (Log)", binwidth=0.2, density_scaling=10)
p.cfdna <- plot.distribution(x="total_cfdna", x_label="Total cfDNA", binwidth=0.05, density_scaling=10)
p.cfdna_log <- plot.distribution(x="total_cfdna_log", x_label="Total cfDNA (Log)", binwidth=0.05, density_scaling = 5)

p.cox3_grid <- grid.arrange(p.cox3, p.cox3_log, nrow=1)
ggsave(paste0(output_dir, "cox3_combined.png"), p.cox3_grid, dpi = 300, width = 4000, height = 1500, units = "px")

p.nd2_grid <- grid.arrange(p.nd2, p.nd2_log, nrow=1)
ggsave(paste0(output_dir, "nd2_combined.png"), p.nd2_grid, dpi = 300, width = 4000, height = 1500, units = "px")

p.cfdna_grid <- grid.arrange(p.cfdna, p.cfdna_log, nrow=1)
ggsave(paste0(output_dir, "total_cfdna_combined.png"), p.cfdna_grid, dpi = 300, width = 4000, height = 1500, units = "px")

p.cox3nd2 <- grid.arrange(p.cox3, p.cox3_log, p.nd2, p.nd2_log, nrow=2)
ggsave(paste0(output_dir, "cox3_nd2_combined.png"), p.cox3nd2, dpi = 300, width = 4000, height = 3000, units = "px")

# Log transformed total_cfdna does not have a normal distribution either
# test_cfdna_distribution <- subset(merged_df, total_cfdna_log != -Inf)
# shapiro.test(test_cfdna_distribution$total_cfdna_log)

# COX3 against ND2 Correlation
ggpubr::ggscatter(merged_df, x = "nd2_log", y = "cox3_log", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = FALSE, cor.method = "spearman", cor.coef.size=8) +
  xlab("Log ND2") +
  ylab("Log COX3") +
  labs(title = "COX3 against ND2") +
  theme_options
ggsave(paste0(output_dir, "cox3_nd2_correlation.png"), dpi = 300, width = 2500, height = 2000, units = "px")

ggpubr::ggscatter(merged_df, x = "total_cfdna_log", y = "cox3_log", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = FALSE, cor.method = "spearman", cor.coef.size=8) +
  xlab("Log total cfDNA") +
  ylab("Log COX3") +
  labs(title = "COX3 against total cfDNA") +
  theme_options
ggsave(paste0(output_dir, "cox3_total_cfdna_correlation.png"), dpi = 300, width = 2500, height = 2000, units = "px")

print("Script successfully completed.")
