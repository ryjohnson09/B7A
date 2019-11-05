# heatmaps.R
# Ryan Johnson
# 5 Nov 2019
# Produce relative abundance heatmaps from 
#  B7A microbiome data

## Load Packages -------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(ggthemes)
library(RColorBrewer)



## Load Data and set Variables -----------------------------------------------
# Data load
ps2 <- readRDS("results/microbiome_analysis/ps2.rds")

# Variables
num_top_ASV <- 30 # How many most abundant ASVs to consider in plots
num_top_genus <- 30 # How many most abundant genera to consider in plots



## Covert to Relative Abundance ----------------------------------------------
# ASV
# Get top N ASVs
topN_ASV <- names(sort(taxa_sums(ps2), decreasing = TRUE))[1:num_top_ASV]

# Convert to relative abundance
ps_ASV_relabun <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU) * 100)

# Extract the top N ASVs
ps_ASV_topN <- prune_taxa(topN_ASV, ps_ASV_relabun)

# Convert ASV abundance to data frame
ASV_abundance_table_topN <- psmelt(ps_ASV_topN)


# Genus
# Combine data to Genus Level
ps_genus <- tax_glom(ps2, "Genus", NArm = TRUE)

# Convert genus data to relative abundance
ps_genus_relabun <- transform_sample_counts(ps_genus, 
                                            function(OTU) OTU/sum(OTU) * 100)

# Get top N genera
topN_genus <- names(sort(taxa_sums(ps_genus), decreasing=TRUE))[1:num_top_genus]

# Extract the top N genera
ps_genus_topN <- prune_taxa(topN_genus, ps_genus_relabun)

# Convert top genus relative abundance to data frame
taxa_abundance_table_genus <- psmelt(ps_genus_topN)



## Heatmap Plots -------------------------------------------------------------------
# Get order of most to least abundance top ASVs
top_ASV_N <- ASV_abundance_table_topN %>% 
  group_by(OTU, Genus) %>% 
  summarise(mean_abun = mean(Abundance)) %>% 
  arrange(desc(mean_abun))

# Refactor
ASV_abundance_table_topN$OTU <- factor(ASV_abundance_table_topN$OTU,
                                       levels = rev(top_ASV_N$OTU))

###################################################################

# Heatmap of ASVs mod sev disease
ASV_heatmap_mod_sev <- ASV_abundance_table_topN %>% 
  ggplot(aes(fill = Abundance, y = OTU, x = Study_ID)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c(breaks = c(0, 0.002, 0.02, 0.2, 2, 20), 
                       labels = c(0, 0.002, 0.02, 0.2, 2, 20),
                       trans = scales::pseudo_log_trans(sigma = 0.002),
                       limits = c(0,20)) +
  scale_y_discrete(limits = as.character(top_ASV_N$OTU), 
                   labels = as.character(top_ASV_N$Genus)) +
  labs(x = "",
       y = "Genus of ASV",
       fill = "Relative Abundance",
       title = paste0("ASV (top ", num_top_ASV, ") Relative Abundance"),
       subtitle = "Ordered by Day and Disease Severity",
       caption = "Y: Mod-Severe Disease\nN: Mild Disease") +
  facet_wrap(~Study_Day + Mod.Sev.Diarrhea, scales = "free_x", nrow = 1) +
  theme(
    legend.key.height = unit(0.5, "in"),
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, "in")
  )

ggsave(filename = "results/microbiome_analysis/heatmaps/ASV_heatmap_mod_sev.png",
       plot = ASV_heatmap_mod_sev,
       height = 8,
       width = 18)

###################################################################

# Get order of most to least abundance top genera
top_genera_N <- taxa_abundance_table_genus %>% 
  group_by(Genus) %>% 
  summarise(mean_abun = mean(Abundance)) %>% 
  arrange(desc(mean_abun))

# Refactor
taxa_abundance_table_genus$Genus <- factor(taxa_abundance_table_genus$Genus,
                                             levels = rev(top_genera_N$Genus))

###################################################################

# Genus heatmap mod sev
genus_heatmap_mod_sev <- taxa_abundance_table_genus %>% 
  ggplot(aes(fill = Abundance, y = Genus, x = Study_ID)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c(breaks = c(0, 0.01, 0.1, 1, 10, 100), 
                       labels = c(0, 0.01, 0.1, 1, 10, 100),
                       trans = scales::pseudo_log_trans(sigma = 0.001),
                       limits = c(0,100)) +
  labs(x = "",
       y = "Genus",
       fill = "Relative Abundance",
       title = paste0("Genus (top ", num_top_ASV, ") Relative Abundance"),
       subtitle = "Ordered by Day and Disease Severity",
       caption = "Y: Mod-Severe Disease\nN: Mild Disease") +
  facet_wrap(~Study_Day + Mod.Sev.Diarrhea, scales = "free_x", nrow = 1) +
  theme(
    legend.key.height = unit(0.5, "in"),
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, "in")
  )

ggsave(filename = "results/microbiome_analysis/heatmaps/genus_heatmap_mod_sev.png",
       plot = genus_heatmap_mod_sev,
       height = 8,
       width = 18)

#####################################################################

# Genus heatmap mod sev mean values
genus_heatmap_mean_mod_sev <- taxa_abundance_table_genus %>% 
  ungroup() %>% 
  group_by(Genus, Study_Day, Mod.Sev.Diarrhea) %>% 
  summarize(mean_abundance = mean(Abundance)) %>% 
  ggplot(aes(fill = mean_abundance, y = Genus, x = Mod.Sev.Diarrhea)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c(breaks = c(0, 0.01, 0.1, 1, 10, 100), 
                       labels = c(0, 0.01, 0.1, 1, 10, 100),
                       trans = scales::pseudo_log_trans(sigma = 0.001),
                       limits = c(0,100),
                       option = "B") +
  labs(y = "Genus",
       x = "Moderate to Severe Diarrhea?",
       fill = "Relative Abundance",
       title = paste0("Genus (top ", num_top_ASV, ") Mean Relative Abundance")) +
  facet_wrap(~Study_Day, scales = "free_x", nrow = 1) +
  theme(
    axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5)
  )

ggsave(filename = "results/microbiome_analysis/heatmaps/genus_heatmap_mean_mod_sev.png",
       plot = genus_heatmap_mean_mod_sev,
       height = 8,
       width = 10)

#####################################################################

# Genus heatmap mean naive
genus_heatmap_mean_naive <- taxa_abundance_table_genus %>% 
  ungroup() %>% 
  group_by(Genus, Study_Day, Naïve) %>% 
  summarize(mean_abundance = mean(Abundance)) %>% 
  ggplot(aes(fill = mean_abundance, y = Genus, x = Naïve)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c(breaks = c(0, 0.01, 0.1, 1, 10, 100), 
                       labels = c(0, 0.01, 0.1, 1, 10, 100),
                       trans = scales::pseudo_log_trans(sigma = 0.001),
                       limits = c(0,100),
                       option = "B") +
  labs(y = "Genus",
       x = "Naive?",
       fill = "Relative Abundance",
       title = paste0("Genus (top ", num_top_ASV, ") Mean Relative Abundance")) +
  facet_wrap(~Study_Day, scales = "free_x", nrow = 1) +
  theme(
    axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5)
  )

ggsave(filename = "results/microbiome_analysis/heatmaps/genus_heatmap_mean_naive.png",
       plot = genus_heatmap_mean_naive,
       height = 8,
       width = 10)

#####################################################################

# Genus heatmap pre post abx
taxa_abundance_table_genus$pre_post_abx <- factor(
  taxa_abundance_table_genus$pre_post_abx, levels = c("Pre", "Post")
)

genus_heatmap_mean_abx <- taxa_abundance_table_genus %>% 
  ungroup() %>% 
  group_by(Genus, Study_Day, pre_post_abx) %>% 
  summarize(mean_abundance = mean(Abundance)) %>% 
  ggplot(aes(fill = mean_abundance, y = Genus, x = pre_post_abx)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c(breaks = c(0, 0.01, 0.1, 1, 10, 100), 
                       labels = c(0, 0.01, 0.1, 1, 10, 100),
                       trans = scales::pseudo_log_trans(sigma = 0.001),
                       limits = c(0,100),
                       option = "B") +
  labs(y = "Genus",
       x = "Pre or Post Antibiotics",
       fill = "Relative Abundance",
       title = paste0("Genus (top ", num_top_ASV, ") Mean Relative Abundance")) +
  facet_wrap(~Study_Day, scales = "free_x", nrow = 1) +
  theme(
    axis.text.x = element_text(angle = 0, size = 8, vjust = 0.5)
  )

ggsave(filename = "results/microbiome_analysis/heatmaps/genus_heatmap_mean_abx.png",
       plot = genus_heatmap_mean_abx,
       height = 8,
       width = 10)