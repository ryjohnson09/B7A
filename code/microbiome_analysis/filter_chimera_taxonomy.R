# filter_chimera_taxonomy.R
# Ryan Johnson
# 4 November 2019
# Remove improper length reads, remove chimeras 
#  and assign taxonomy from the B7A microbiome 
#  data using DADA2

# Load Libraries --------------------------------------------------
library(dada2)
library(tidyverse)


## Load in DADA2 data ----------------------------------------------
seqtab <- readRDS("results/microbiome_analysis/seqtab.rds")


## Remove sequences that are incorrect length ---------------------------------
# Based on the primers used, I'm expecting fragments in the 465 bp size range,
# Filter everything below 438 and 
#  above 472 (based on `table(nchar(getSequences(seqtab)))`)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 438:472]


## Remove Chimeras ------------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab2, 
                                    method="consensus", 
                                    multithread=TRUE, 
                                    verbose = TRUE)

chimera_removed_dimensions <- dim(seqtab.nochim)

chimera_removed_divided_by_original <- sum(seqtab.nochim)/sum(seqtab2)

# Save the results to text file
write_lines(x = chimera_removed_divided_by_original, 
           path = "results/microbiome_analysis/chimera_remaining.txt")


## Assign Taxonomy ------------------------------------------------------------

tax <- assignTaxonomy(seqs = seqtab.nochim, 
                      refFasta = "data/processed/microbiome_analysis/silva_nr_v132_train_set.fa.gz", 
                      multithread=TRUE)

## Save Results to disk --------------------------------------------------------
saveRDS(seqtab.nochim, "results/microbiome_analysis/seqtab_final.rds")
saveRDS(tax, "results/microbiome_analysis/tax_final.rds")