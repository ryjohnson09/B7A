# abundance_barplots.R
# Ryan Johnson
# 5 Nov 2019
# Analyze polished microbiome data to produce 
#  relative abundance plots at phylum and genus
#  level

## Load Packages -------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(ggthemes)
library(RColorBrewer)



## Load Data and set Variables -----------------------------------------------
# Data load
ps2 <- readRDS("results/microbiome_analysis/ps2.rds")

# Variables
num_top_genus <- 12 # How many most abundant genera to consider in plots



## Covert to Relative Abundance ----------------------------------------------
# Combine data to Phylum Level
ps_phylum <- tax_glom(ps2, "Phylum", NArm = TRUE)

# Convert phylum data to Relative Abundacne
ps_phylum_relabun <- transform_sample_counts(ps_phylum, 
                                             function(OTU) OTU/sum(OTU) * 100)
# Convert phylum relative abundance to data frame
taxa_abundance_table_phylum <- psmelt(ps_phylum_relabun)

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



## Phylum Level Bar Plots -----------------------------------------------------
# Phylum Level By Patient
barplot_phylum_patient <- taxa_abundance_table_phylum %>% 
  mutate(Study_Day = factor(Study_Day)) %>% 
  ggplot(aes(x = Study_Day, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Study Day",
       y = "Relative Abundance",
       title = "Relative Abundance at Phylum Level") +
  facet_wrap(~ Study_ID) +
  scale_fill_colorblind() +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = "results/microbiome_analysis/barplots/barplot_phylum_patient.png", 
       plot = barplot_phylum_patient, 
       height = 8,
       width = 15)

#############################################################

# Phylum abundance by time point
# Get order of study_id based on Proteobacteria abundance
study_ID_order_proteo <- taxa_abundance_table_phylum %>% 
  filter(Phylum == "Proteobacteria") %>% 
  arrange(desc(Abundance)) %>% 
  pull(Sample)

barplot_phylum_time <- taxa_abundance_table_phylum %>% 
  mutate(Study_Day = factor(Study_Day)) %>% 
  ggplot(aes(x = factor(Sample, levels = study_ID_order_proteo), y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "",
       y = "Relative Abundance",
       title = "Relative Abundance at Phylum Level",
       subtitle = "Samples ordered by Proteobacteria Abundance") +
  facet_wrap(~ Study_Day, scales = "free") +
  scale_fill_colorblind() +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, units = "in"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = "results/microbiome_analysis/barplots/barplot_phylum_time.png", 
       plot = barplot_phylum_time, 
       height = 8,
       width = 15)

#############################################################

# Add in Moderate to Severe Disease
taxa_abundance_table_phylum$Mod.Sev.Diarrhea <- 
  factor(taxa_abundance_table_phylum$Mod.Sev.Diarrhea, ordered = TRUE)

barplot_phylum_Mod_Sev <- taxa_abundance_table_phylum %>% 
  mutate(Study_Day = factor(Study_Day)) %>% 
  ggplot(aes(x = Study_Day, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Study Day",
       y = "Relative Abundance",
       title = "Relative Abundance at Phylum Level",
       caption = "Y: Mod-Severe Disease\nN: Mild Disease") +
  facet_wrap(~ Mod.Sev.Diarrhea + Study_ID) +
  scale_fill_colorblind() +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = "results/microbiome_analysis/barplots/barplot_phylum_mod_sev.png", 
       plot = barplot_phylum_Mod_Sev, 
       height = 10,
       width = 15)

#############################################################

# Phylum Naive vs Non-Naive
barplot_phylum_naive <- taxa_abundance_table_phylum %>% 
  mutate(Study_Day = factor(Study_Day)) %>% 
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "",
       y = "Relative Abundance",
       title = "Relative Abundance at Phylum Level",
       caption = "Y: Naive\nN: Not Naive") +
  facet_wrap(~ Study_Day + Naïve, scales = "free") +
  scale_fill_colorblind() +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, units = "in"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = "results/microbiome_analysis/barplots/barplot_phylum_naive.png", 
       plot = barplot_phylum_naive, 
       height = 10,
       width = 15)

#############################################################

# Phylum pre-post antibiotics
taxa_abundance_table_phylum$pre_post_abx <- factor(
  taxa_abundance_table_phylum$pre_post_abx, levels = c("Pre", "Post")
)

barplot_phylum_abx <- taxa_abundance_table_phylum %>% 
  mutate(Study_Day = factor(Study_Day)) %>% 
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "",
       y = "Relative Abundance",
       title = "Relative Abundance at Phylum Level",
       caption = "") +
  facet_wrap(~ Study_Day + pre_post_abx, scales = "free") +
  scale_fill_colorblind() +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, units = "in"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = "results/microbiome_analysis/barplots/barplot_phylum_abx.png", 
       plot = barplot_phylum_abx, 
       height = 10,
       width = 15)



## Genus Level Bar Plots -----------------------------------------------------
# Get fill colors
getPalette <- colorRampPalette(brewer.pal(num_top_genus, "Paired"))

# Genus level by patient
barplot_genus_patient <- taxa_abundance_table_genus %>% 
  mutate(Study_Day = factor(Study_Day)) %>% 
  ggplot(aes(x = Study_Day, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Study Day",
       y = "Relative Abundance",
       title = paste0("Genus (top ", num_top_genus, ") Relative Abundance")) +
  facet_wrap(~ Study_ID) +
  scale_fill_manual(values = getPalette(num_top_genus)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = "results/microbiome_analysis/barplots/barplot_genus_patient.png", 
       plot = barplot_genus_patient, 
       height = 8,
       width = 15)

#############################################################

# Genus by patient and mod severe disease
barplot_genus_mod_sev <- taxa_abundance_table_genus %>% 
  mutate(Study_Day = factor(Study_Day)) %>% 
  ggplot(aes(x = Study_Day, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Study Day",
       y = "Relative Abundance",
       title = paste0("Genus (top ", num_top_genus, ") Relative Abundance"),
       caption = "Y: Mod-Severe Disease\nN: Mild Disease") +
  facet_wrap(~ Mod.Sev.Diarrhea + Study_ID) +
  scale_fill_manual(values = getPalette(num_top_genus)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = "results/microbiome_analysis/barplots/barplot_genus_mod_sev.png", 
       plot = barplot_genus_mod_sev, 
       height = 10,
       width = 15)

#############################################################

# Genus abundance by day
## Get order of samples based on Escherichia abundance
Sample_ID_order_ecoli <- taxa_abundance_table_genus %>% 
  filter(Genus == "Escherichia/Shigella") %>% 
  arrange(desc(Abundance)) %>% 
  pull(Sample)

barplot_genus_time <- taxa_abundance_table_genus %>% 
  mutate(Study_Day = factor(Study_Day)) %>% 
  ggplot(aes(x = factor(Sample, levels = Sample_ID_order_ecoli), 
             y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "",
       y = "Relative Abundance",
       title = paste0("Genus (top ", num_top_genus, ") Relative Abundance"),
       subtitle = "Samples ordered by Escherichia/Shigella Abundance") +
  facet_wrap(~Study_Day, scales = "free") +
  scale_fill_manual(values = getPalette(num_top_genus)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, units = "in"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = "results/microbiome_analysis/barplots/barplot_genus_time.png", 
       plot = barplot_genus_time, 
       height = 10,
       width = 15)

#############################################################

# Genus Naive vs Non-Naive
barplot_genus_naive <- taxa_abundance_table_genus %>% 
  mutate(Study_Day = factor(Study_Day)) %>% 
  ggplot(aes(x = factor(Sample, levels = Sample_ID_order_ecoli),
             y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "",
       y = "Relative Abundance",
       title = paste0("Genus (top ", num_top_genus, ") Relative Abundance"),
       caption = "Y: Naive\nN: Not Naive") +
  facet_wrap(~ Study_Day + Naïve, scales = "free") +
  scale_fill_manual(values = getPalette(num_top_genus)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, units = "in"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = "results/microbiome_analysis/barplots/barplot_genus_naive.png", 
       plot = barplot_genus_naive, 
       height = 10,
       width = 15)

#############################################################

# Genus pre-post antibiotics
taxa_abundance_table_genus$pre_post_abx <- factor(
  taxa_abundance_table_genus$pre_post_abx, levels = c("Pre", "Post")
)

barplot_genus_abx <- taxa_abundance_table_genus %>% 
  mutate(Study_Day = factor(Study_Day)) %>% 
  ggplot(aes(x = factor(Sample, levels = Sample_ID_order_ecoli),
             y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "",
       y = "Relative Abundance",
       title = paste0("Genus (top ", num_top_genus, ") Relative Abundance"),
       caption = "") +
  facet_wrap(~ Study_Day + pre_post_abx, scales = "free") +
  scale_fill_manual(values = getPalette(num_top_genus)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, units = "in"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = "results/microbiome_analysis/barplots/barplot_genus_abx.png", 
       plot = barplot_genus_abx, 
       height = 10,
       width = 15)
