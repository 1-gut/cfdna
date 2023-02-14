# Phyloseq Downstream Analysis from processed files from cfdna_nextflow
# ------------------------------------------------------------------------

output_dir <- "output/phyloseq/"

library(dplyr)
library(ggplot2)
library(phyloseq)
library(stringr)
library(vegan)
library(DESeq2)
library(scales)
library(RColorBrewer)
library(ggsignif)
library(gridExtra)
source("./theme_options.R")
# ------------------------------------------------------------------------
# SWITCH DATA SOURCES BETWEEN DIRECT AND UNMAPPED
phyloseq_metadata <- read.csv("data/phyloseq_metadata.csv", fileEncoding = "UTF-8", row.names=4)

phyloseq_metadata$ibd_status_collapsed <- as.factor(phyloseq_metadata$ibd_status_collapsed)
phyloseq_metadata$ibd_status_collapsed <- relevel(phyloseq_metadata$ibd_status_collapsed, "Remission")
phyloseq_metadata$ibd_status_collapsed <- relevel(phyloseq_metadata$ibd_status_collapsed, "Active")

#SWITCH between running on unmapped sequences or raw sequences
biom_data <- import_biom("data/03052022/reports_pipeline_v1.3/kraken2_unmapped/biom/collated_kraken_bracken_unmapped.biom")
# biom_data <- import_biom("data/03052022/reports_pipeline_v1.3/kraken2/biom/collated_kraken_bracken.biom")

# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Initial pre-processing
# ------------------------------------------------------------------------
# Merge biom file with annotated metadata
sampledata <- sample_data(phyloseq_metadata)
biom_data = merge_phyloseq(biom_data, sampledata)

# Replace raw Kraken2 taxonomic ranks with descriptive ranks
biom_data@tax_table@.Data <- substring(biom_data@tax_table@.Data, 4)
colnames(biom_data@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Rename samples to just the ID.
sample_names(biom_data) <- sample_data(biom_data)$study_id


# ------------------------------------------------------------------------
# Plot raw reads from biom files
# ------------------------------------------------------------------------

# Show number of reads in each sample and plot
sample_sums(biom_data)
deep <- data.frame(Samples = sample_names(biom_data),
                   Reads = sample_sums(biom_data),
                   Grouping = sample_data(biom_data)$ibd_status_collapsed)
ggplot(data = deep, mapping = aes(x = Samples,y = Reads)) +
  geom_col() +
  scale_y_continuous(label = comma) +
  theme_options

deep %>%
  subset(Samples != "PC") %>%
  ggplot(aes(x=Samples, y=Reads)) +
  geom_col() +
  scale_y_continuous(label = comma) +
  labs(title="Unmapped reads (after removing host reads)") +
  theme_options
ggsave(paste0(output_dir, "unmapped_reads.png"),
       dpi = 300, width = 3000, height = 2000, units = "px")


# ------------------------------------------------------------------------
# Compute Rarefaction Curves
# ------------------------------------------------------------------------

# Computes rarefaction curve with and without PC

# First, compute the rarefaction values (not informative)
# rare_results_with_pc <- vegan::rarecurve(t(otu_table(biom_data)), step=5000, cex=0.8, tidy=TRUE)

# Remove PC and compute values
rarefaction_without_pc <- subset_samples(biom_data, ibd_status_collapsed != 'Positive control')
rare_results <- vegan::rarecurve(t(otu_table(rarefaction_without_pc)), step=5000, cex=0.8, tidy=TRUE)

# Plot rarefaction curves with ggplot
rare_results %>%
  ggplot(aes(x=Sample, y=Species, group=Site, col=Site)) +
  geom_line() +
  geom_label(aes(label = Site),
             data = rare_results %>% group_by(Site) %>% filter(Sample == max(Sample)),
             nudge_x = 0.35,
             size = 4) +
  labs(title="Rarefaction Curve") +
  scale_x_continuous(label = comma) +
  theme_options
ggsave(paste0(output_dir, "rarefaction_curve.png"), 
       dpi = 300, width = 2000, height = 1500, units = "px")

# Rarefaction curve with PC is not informative
# rare_results_with_pc %>%
#   ggplot(aes(x=Sample, y=Species, group=Site, col=Site)) +
#   geom_line() +
#   geom_label(aes(label = Site),
#              data = rare_results_with_pc %>% group_by(Site) %>% filter(Sample == max(Sample)),
#              nudge_x = 0.35,
#              size = 4) +
#   labs(title="Rarefaction Curve with PC") +
#   scale_x_continuous(label = comma) +
#   theme_options
# ggsave("output/phyloseq/rarefaction_curve_with_pc.png", 
#        dpi = 300, width = 2000, height = 1500, units = "px")


# ------------------------------------------------------------------------
# Removes data with empty Genus before percentage transformation
# ------------------------------------------------------------------------

# Remove empty genus in preparation for abundance visualization
summary(biom_data@tax_table@.Data== "")
biom_data <- subset_taxa(biom_data, Genus != "")

# add percentages calculation for beta diversity plotting to normalise reads across samples
percentages  = transform_sample_counts(biom_data, function(x) x*100 / sum(x) )


# ------------------------------------------------------------------------
# Beta Diversity
# ------------------------------------------------------------------------
meta.ord <- ordinate(physeq = percentages, method = "NMDS", 
                     distance = "bray")

# Beta Diversity by Grouping
# plot_ordination(physeq = percentages, ordination = meta.ord, color="study_group_name") +
#   geom_point(size=3) +
#   labs(title="Beta Diversity by Grouping") +
#   theme_options +
#   theme(legend.position="right")

# Beta diversity by IBD Activity
plot_ordination(physeq = percentages, ordination = meta.ord, color="ibd_status_collapsed") +
  labs(
    title="Beta Diversity by IBD Activity",
    caption="Distance calculated using Bray-Curtis method",
    color="IBD Status") +
  geom_point(size=4) +
  theme_options +
  theme(legend.position="right") +
  scale_fill_brewer(palette="Pastel1")
ggsave(
  paste0(output_dir, "beta_diversity_by_ibd_status_collapsed.png"),
  width = 3000, height = 2000, units = "px", limitsize = FALSE
)


# ANOSIM Test


otu_table <- t(as.data.frame(otu_table(biom_data)))
otu_table_cleaned <- data.frame(study_id = row.names(otu_table), otu_table)
otu_table_metadata <- phyloseq_metadata %>% select(c("study_id", "ibd_status_collapsed_code"))
otu_table_cleaned <- left_join(otu_table_cleaned, otu_table_metadata, by="study_id")
ano = anosim(otu_table, otu_table_cleaned$ibd_status_collapsed_code, distance = "bray", permutations = 9999)
ano

# Permanova test

permanova = adonis2(otu_table ~ ibd_status_collapsed_code, data=otu_table_metadata, permutations=999, method="bray")
permanova
# Plot beta diversity by activity with split plot of phylum
# plot_ordination(percentages, meta.ord, type="split", color="Phylum", shape="ibd_status_collapsed", title="Phylum") + theme_options

# ------------------------------------------------------------------------
# Raw and relative abundance visualization
# ------------------------------------------------------------------------
# Raw abundance and relative abundance plots per sample by Phylum


plot_abundance <- function(taxonomic_rank, abundance_threshold=0.01) {
  # Aggregate species at the specified taxonomic_rank for both raw and percentage datasets
  glom <- phyloseq::tax_glom(percentages, taxrank = taxonomic_rank)
  percentages_data <- phyloseq::psmelt(glom)
  
  # raw reads not really informative
  # raw <- tax_glom(physeq = biom_data, taxrank = taxonomic_rank)
  # raw_data <- phyloseq::psmelt(raw)
  # raw_plot <- raw_data %>%
  #   subset(OTU!=9606) %>% # remove human reads
  #   subset(Abundance > 30) %>% # min 10 raw reads
  #   ggplot(aes(x=Sample, y=Abundance, fill=.data[[taxonomic_rank]]))+ 
  #   geom_bar(aes(), stat="identity", position="stack", show.legend=TRUE) +
  #   labs(title=paste0("Raw Reads by ", taxonomic_rank)) +
  #   theme_options +
  #   theme(legend.position = "right")
  # 
  # print(raw_plot)

  rel_plot <- percentages_data %>%
    subset(OTU!=9606) %>%
    subset(Abundance > abundance_threshold) %>%
    ggplot(aes(x=Sample, y=Abundance, fill=.data[[taxonomic_rank]]))+ 
    geom_bar(aes(), stat="identity", position="stack", show.legend=TRUE) +
    labs(title=paste0("Relative Reads by ", taxonomic_rank)) +
    theme_options +
    theme(legend.position = "right")
  
  print(rel_plot)
  ggsave(paste0(output_dir, "relative_reads_by_", taxonomic_rank, ".png"),
         dpi=300, width=4000, height=2500, units="px")
}

plot_abundance("Phylum")
# plot_abundance("Class")
# plot_abundance("Order", abundance_threshold=0.5)
# plot_abundance("Family", abundance_threshold=0.5)
plot_abundance("Genus", abundance_threshold=0.15)

# ------------------------------------------------------------------------
# Alpha diversity plots
# ------------------------------------------------------------------------

selected_measures <- c("Chao1", "Shannon") # Observed looks similar to Chao1

# Individual Sample Visualization (not informative)
# plot_richness(biom_data, measures=selected_measures) +
#   labs(title="Alpha diversity for each sample", x="Samples") +
#   theme_options
# ggsave(
#   "output/phyloseq/alpha_diversity_by_sample.png",
#   dpi=300, width = 2000, height = 1500, units = "px", limitsize = FALSE
# )

# p <- plot_richness(biom_data, measures=selected_measures, x="study_group_name") +
#   xlab("Group") +
#   labs(title="Alpha diversity by group")
# p$layers <- p$layers[-1]
# p + geom_boxplot() + geom_jitter()
# ggsave(
#   "output/phyloseq/alpha_diversity_by_group.png",
#   dpi=300, width = 2000, height = 1500, units = "px", limitsize = FALSE
# )

# Includes controls
p <- plot_richness(biom_data, measures=selected_measures, x="ibd_status_collapsed") +
  labs(title="Reduced alpha diversity in active disease", x="IBD Status")
p$layers <- p$layers[-1]
p$layers <- p$layers[-1]
p + stat_boxplot(geom="errorbar", width=0.2) + geom_boxplot(width=0.3) + theme_options
# ggsave(
#   "output/phyloseq/alpha_diversity_by_ibd_status_collapsed.png",
#   dpi=300, width = 2000, height = 1500, units = "px", limitsize = FALSE
# )

# Alpha Diversity Chao1 and Shannon Final Plots
remove_controls <- subset_samples(biom_data, ibd_status_collapsed != 'Negative control')
remove_controls <- subset_samples(remove_controls, ibd_status_collapsed != 'Positive control')
p <- plot_richness(remove_controls, measures="Chao1", x="ibd_status_collapsed") +
  xlab("IBD Status")
p$layers <- p$layers[-1]
p$layers <- p$layers[-1]
p_chao1 <- p + stat_boxplot(geom="errorbar", width=0.2) +
  geom_boxplot(width=0.3) +
  theme_options +
  geom_signif(comparisons = list(c("Active", "Healthy control")),
              y_position = 5600,
              annotations = "p = 0.07",
              tip_length = 0.01,
              vjust = 0) +theme(title = element_blank())
p_chao1

p <- plot_richness(remove_controls, measures=c("Shannon"), x="ibd_status_collapsed") +
  xlab("IBD Status")
p$layers <- p$layers[-1]
p$layers <- p$layers[-1]
p_shannon <- p + stat_boxplot(geom="errorbar", width=0.2) +
  geom_boxplot(width=0.3) +
  theme_options +
  geom_signif(comparisons = list(c("Active", "Healthy control")),
              y_position = 0.9,
              annotations = "p = 0.43",
              tip_length = 0.01,
              vjust = 0) +theme(title = element_blank())

p_shannon



p_alpha_diversity_final <- gridExtra::grid.arrange(
  p_chao1, p_shannon, nrow=1
)

ggsave(
  paste0(output_dir, "alpha_diversity_by_ibd_activity.png"),
  p_alpha_diversity_final,
  dpi=300, width = 4000, height = 1500, units = "px", limitsize = FALSE
)


# ------------------------------------------------------------------------
# Alpha diversity stats
# ------------------------------------------------------------------------
alpha_diversity <- data.frame(
  "Observed" = phyloseq::estimate_richness(biom_data, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(biom_data, measures = "Shannon"),
  "Chao1" = phyloseq::estimate_richness(biom_data, measures = "Chao1"),
  "ibd_status_collapsed" = phyloseq::sample_data(biom_data)$ibd_status_collapsed)

summary_table <- alpha_diversity %>%
  subset(ibd_status_collapsed != 'Negative control') %>%
  subset(ibd_status_collapsed != 'Positive control') %>%
  group_by(ibd_status_collapsed) %>%
  dplyr::summarise(median_observed = median(Observed),
            median_shannon = median(Shannon),
            median_chao1 = median(Chao1.Chao1))

alpha_diversity_without_controls <- alpha_diversity %>%
  subset(ibd_status_collapsed != 'Negative control') %>%
  subset(ibd_status_collapsed != 'Positive control')

alpha_diversity.anova.obv <- aov(Observed ~ ibd_status_collapsed, alpha_diversity_without_controls) 
summary(alpha_diversity.anova.obv)
tukey.test.ob <- TukeyHSD(alpha_diversity.anova.obv)
tukey.test.ob

alpha_diversity.anova.chao <- aov(Chao1.Chao1 ~ ibd_status_collapsed, alpha_diversity_without_controls) 
summary(alpha_diversity.anova.chao)
TukeyHSD(alpha_diversity.anova.chao)
# Chao1 anova p=0.07 between HC and active. Post-hoc Rest all non significant

alpha_diversity.anova.shannon <- aov(Shannon ~ ibd_status_collapsed, alpha_diversity_without_controls) 
summary(alpha_diversity.anova.shannon)
tukey.test.shannon <- TukeyHSD(alpha_diversity.anova.shannon)
tukey.test.shannon

# ------------------------------------------------------------------------
# Plot heatmap of top 30 species
# ------------------------------------------------------------------------
sorted_sample_metadata <- arrange(phyloseq_metadata, ibd_status_collapsed)
sort_order <- sorted_sample_metadata$study_id

heatmap_without_humans <- subset_taxa(biom_data, Genus!="Homo")

heatmap_without_humans_glom <- phyloseq::tax_glom(heatmap_without_humans, taxrank = "Genus")

biom_heatmap <- prune_taxa(names(sort(taxa_sums(heatmap_without_humans_glom),TRUE)[1:30]), heatmap_without_humans_glom)
plot_heatmap(biom_heatmap, "nmds", "bray", sample.label="ibd_status", taxa.label="Genus", sample.order=sort_order)
# ggsave(
#   "output/phyloseq/heatmap_top_30.png",
#   width = 2000, height = 1500, units = "px", limitsize = FALSE
# )


hetamap_without_controls <- subset_taxa(biom_data, Genus != "Homo")
heatmap_without_controls <- subset_samples(heatmap_without_humans_glom, ibd_status != "Positive control")
heatmap_without_controls <- subset_samples(heatmap_without_humans_glom, ibd_status != "Negative control")
heatmap_without_controls_glom <- phyloseq::tax_glom(heatmap_without_controls, taxrank = "Genus") # collapse tax ranks to Genus level
biom_heatmap_without_controls <- prune_taxa(names(sort(taxa_sums(heatmap_without_controls_glom),TRUE)[1:30]), heatmap_without_controls_glom) # take topn 30
plot_heatmap(biom_heatmap_without_controls, "nmds", "bray", sample.label="ibd_status", taxa.label="Genus", sample.order=sort_order)

# ggsave(
#   "output/phyloseq/heatmap_without_pc.png",
#   width = 2000, height = 1500, units = "px", limitsize = FALSE
# )

# Only use this to interrogate particular species of interest
# biom_pb <- subset_taxa(biom_data, Phylum=="Proteobacteria")
# plot_heatmap(biom_pb)


# ------------------------------------------------------------------------
# Differential abundance with DeSeq2
# ------------------------------------------------------------------------
diagdds = phyloseq_to_deseq2(biom_data, ~ ibd_status_collapsed_code)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.1
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(biom_data)[rownames(sigtab), ], "matrix"))

head(sigtab)

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

#genus on x-axis
 # ggplot(sigtab, aes(x=log2FoldChange, y=Genus, color=Phylum)) + geom_point(size=6) + 
 #   theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust=1))

print("Script completed.")
