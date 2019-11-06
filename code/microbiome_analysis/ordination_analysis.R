# ordination_analysis.R
# Ryan Johnson
# 6 Nov 2019
# Analyze samples using various ordination 
#  and plot

## Load Packages -------------------------------------------
library(dada2); packageVersion("dada2")
library(tidyverse)
library(RColorBrewer)
library(phyloseq)
library(plotrix)
library(ggthemes)
library(plotly)
library(ggrepel)
theme_set(theme_bw())

## Import Data ---------------------------------------------------------
ps2 <- readRDS("results/microbiome_analysis/ps2.rds")


## Distance Matrix Calculation -----------------------------------------
# Transform data to proportions as appropriate for Bray-Curtis distances
ps_prop <- transform_sample_counts(ps2, function(otu) otu/sum(otu))

# Perform distance calculation (NMDS, PCoA)
ord_nmds_bray <- ordinate(ps_prop, method="NMDS", distance="bray")
ord_pcoa_bray <- ordinate(ps_prop, method="PCoA", distance="bray")

# Convert to tibbles
ord_nmds_bray_data <- as_tibble(plot_ordination(ps_prop, ord_nmds_bray)$data, rownames = "Sample")
ord_pcoa_bray_data <- as_tibble(plot_ordination(ps_prop, ord_pcoa_bray)$data, rownames = "Sample")

# Factor
ord_nmds_bray_data$Mod.Sev.Diarrhea <- factor(ord_nmds_bray_data$Mod.Sev.Diarrhea, 
                                              levels = c("Y", "N"))
ord_nmds_bray_data$Naïve <- factor(ord_nmds_bray_data$Naïve, 
                                              levels = c("Y", "N"))
ord_nmds_bray_data$pre_post_abx <- factor(ord_nmds_bray_data$pre_post_abx, 
                                   levels = c("Pre", "Post"))
ord_pcoa_bray_data$Mod.Sev.Diarrhea <- factor(ord_nmds_bray_data$Mod.Sev.Diarrhea, 
                                              levels = c("Y", "N"))
ord_pcoa_bray_data$Naïve <- factor(ord_nmds_bray_data$Naïve, 
                                   levels = c("Y", "N"))
ord_pcoa_bray_data$pre_post_abx <- factor(ord_nmds_bray_data$pre_post_abx, 
                                          levels = c("Pre", "Post"))


## Plots -----------------------------------------------------------------
# Get colors
getPalette <- colorRampPalette(brewer.pal(
  length(unique(ord_nmds_bray_data$Study_Day)), "Paired"))

# Bray NMDS
bray_nmds_day <- ggplot(data = ord_nmds_bray_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = factor(Study_Day)), size = 2) +
  labs(title = "Bray-Curtis: NMDS",
       color = "Study Day") +
  scale_color_manual(values = getPalette(length(unique(ord_nmds_bray_data$Study_Day)))) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

ggsave(plot = bray_nmds_day,
       filename = "results/microbiome_analysis/ordination_analysis/bray_nmds_day.png",
       height = 5,
       width = 7)

###################################################

# Bray pcoa
bray_pcoa_day <- ggplot(data = ord_pcoa_bray_data, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = factor(Study_Day)), size = 2) +
  labs(title = "Bray-Curtis: PCoA",
       color = "Study Day") +
  scale_color_manual(values = getPalette(length(unique(ord_nmds_bray_data$Study_Day)))) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

ggsave(plot = bray_pcoa_day,
       filename = "results/microbiome_analysis/ordination_analysis/bray_pcoa_day.png",
       height = 5,
       width = 7)

###################################################

# Bray nmds mod sev 
bray_nmds_mod_sev <- ggplot(data = ord_nmds_bray_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Mod.Sev.Diarrhea), size = 2) +
  labs(title = "Bray-Curtis: NMDS",
       subtitle = "Colored by Mod-Sev Diarrhea",
       color = "Mod to Severe\nDiarrhea?") +
  scale_color_colorblind() +
  facet_wrap(~Study_Day) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

ggsave(plot = bray_nmds_mod_sev,
       filename = "results/microbiome_analysis/ordination_analysis/bray_nmds_mod_sev.png",
       height = 7,
       width = 9)

#######################################################

# Bray pcoa mod sev 
bray_pcoa_mod_sev <- ggplot(data = ord_pcoa_bray_data, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = Mod.Sev.Diarrhea), size = 2) +
  labs(title = "Bray-Curtis: PCoA",
       subtitle = "Colored by Mod-Sev Diarrhea",
       color = "Mod to Severe\nDiarrhea?") +
  scale_color_colorblind() +
  facet_wrap(~Study_Day) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

ggsave(plot = bray_pcoa_mod_sev,
       filename = "results/microbiome_analysis/ordination_analysis/bray_pcoa_mod_sev.png",
       height = 7,
       width = 9)

######################################################

# Bray nmds naive
bray_nmds_naive <- ggplot(data = ord_nmds_bray_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Naïve), size = 2) +
  labs(title = "Bray-Curtis: NMDS",
       subtitle = "Colored by Naive (Yes or No)",
       color = "Naïve?") +
  scale_color_colorblind() +
  facet_wrap(~Study_Day) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

ggsave(plot = bray_nmds_naive,
       filename = "results/microbiome_analysis/ordination_analysis/bray_nmds_naive.png",
       height = 7,
       width = 9)

#######################################################

# Bray pcoa naive
bray_pcoa_naive <- ggplot(data = ord_pcoa_bray_data, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = Naïve), size = 2) +
  labs(title = "Bray-Curtis: PCoA",
       subtitle = "Colored by Mod-Sev Diarrhea",
       color = "Naïve?") +
  scale_color_colorblind() +
  facet_wrap(~Study_Day) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

ggsave(plot = bray_pcoa_naive,
       filename = "results/microbiome_analysis/ordination_analysis/bray_pcoa_naive.png",
       height = 7,
       width = 9)

########################################################

# Bray nmds abx
bray_nmds_abx <- ggplot(data = ord_nmds_bray_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = pre_post_abx), size = 2) +
  labs(title = "Bray-Curtis: NMDS",
       subtitle = "Pre or Post Abx Treatment",
       color = "Pre or Post\nAntibiotics") +
  scale_color_colorblind() +
  facet_wrap(~Study_Day) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

ggsave(plot = bray_nmds_abx,
       filename = "results/microbiome_analysis/ordination_analysis/bray_nmds_abx.png",
       height = 7,
       width = 9)

########################################################

# Bray pcoa abx
bray_pcoa_abx <- ggplot(data = ord_pcoa_bray_data, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = pre_post_abx), size = 2) +
  labs(title = "Bray-Curtis: PCoA",
       subtitle = "Pre or Post Abx Treatment",
       color = "Pre or Post\nAntibiotics") +
  scale_color_colorblind() +
  facet_wrap(~Study_Day) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

ggsave(plot = bray_pcoa_abx,
       filename = "results/microbiome_analysis/ordination_analysis/bray_pcoa_abx.png",
       height = 7,
       width = 9)
