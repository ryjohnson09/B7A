# phyloseq_filtering.R
# Ryan Johnson
# 5 Nov 2019
# Using Phyloseq package to explore and filter 
#  microbiome data

## Load Packages --------------------------------------------------
library(tidyverse)
library(phyloseq)



## Load Data from DADA2 -------------------------------------------
seqtab.nochim <- readRDS("results/microbiome_analysis/seqtab_final.rds")
tax <- readRDS("results/microbiome_analysis/tax_final.rds")



## Create Metadata table ------------------------------------------
# Read in metadata
metadata <- read_csv("data/processed/microbiome_analysis/B7A_metadata_clean.csv") %>% 
  as.data.frame()

# Read in sample names from seqtab_final.rds
samples_out <- rownames(seqtab.nochim)

# Check for missing samples
missing_samples <- setdiff(samples_out, str_replace(metadata$TUBE, "_S", "-S"))

# Give Rownames to metadata
rownames(metadata) <- str_replace(metadata$TUBE, "_S", "-S")

# Add the metadata to the seqtab.nochim data
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(tax))
ps



## Remove unwanted phyla ---------------------------------------------
# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Remove unwanted phyla
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized", "Euglenozoa"))



## Filter by prevalence ---------------------------------------------
# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(ps),
                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(ps),
                     tax_table(ps))

# Compute the total and average prevalences of each feature:
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Subset to the remaining phyla
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))

# Plot total abundance (x axis) vs prevelance among samples (y axis)
prev_plot <-
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps), color = Phylum)) +
  # Include a guess for parameter
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap( ~ Phylum) + theme(legend.position = "none")

ggsave(filename = "results/microbiome_analysis/prev_plot_before_filt.png", plot = prev_plot)

# Define phyla to filter
filterPhyla <- c("Cyanobacteria", "Epsilonbacteraeota", "Fusobacteria", "Lentisphaerae", 
                 "Patescibacteria", "Synergistetes", "Tenericutes", "uncharacterized")

# Filter phyla including unidentified phylum
ps1 <- subset_taxa(ps, !Phylum %in% filterPhyla)

# Define prevalence threshold as 5% of total samples
prevalenceThreshold <- 0.05 * nsamples(ps)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 <- prune_taxa(keepTaxa, ps1)



## Additional Filtering ------------------------------------------------------------
# Remove Sample B7A062 for the time being (since only have -1 sample)
ps2 <- subset_samples(physeq = ps2, Study_ID != "B7A062")



## Write to file --------------------------------------------
saveRDS(ps2, "results/microbiome_analysis/ps2.rds")

