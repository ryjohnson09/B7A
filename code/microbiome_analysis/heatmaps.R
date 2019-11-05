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
