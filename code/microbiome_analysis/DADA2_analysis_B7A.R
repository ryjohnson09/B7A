# DADA2_analysis_B7A
# Ryan Johnson
# 4 November 2019
# Process the B7A 16S rRNA using DADA2

## Install Packages -----------------------------------------------------------
library(dada2); packageVersion("dada2")
suppressMessages(suppressWarnings(library(tidyverse)))
library(phyloseq)
library(tidyverse)


## Set paths to reads and filtered reads--------------------------------------
# Path to Forward and Reverse fastq files
pathF <- "data/raw/microbiome_sequencing/ETEC/ETEC_sequences/F_reads/" 
pathR <- "data/raw/microbiome_sequencing/ETEC/ETEC_sequences/R_reads/"

# Where to place filtered files
# Creat file if not already present
if (!dir.exists("data/processed/microbiome_analysis/filtered_F_reads")){
  dir.create("data/processed/microbiome_analysis/filtered_F_reads")
  dir.create("data/processed/microbiome_analysis/filtered_R_reads")
}

# Set path to filtered reads
filtpathF <- "data/processed/microbiome_analysis/filtered_F_reads"
filtpathR <- "data/processed/microbiome_analysis/filtered_R_reads"

# Get list of fastq files
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))

# Extract sample names
sample_names <- paste(
  sapply(strsplit(basename(fastqFs), "_"), `[`, 1), 
  sapply(strsplit(basename(fastqFs), "_"), `[`, 2), 
  sep = "-")

# Ensure equal number of forward and reverse reads
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")


## Check read quality ------------------------------------------------------
# Extract subset of reads from F and R and assess read quality.
# Get full path of raw reads (F and R)
fastqFs_path <- sort(list.files(pathF, pattern="fastq.gz", full.names = TRUE))
fastqRs_path <- sort(list.files(pathR, pattern="fastq.gz", full.names = TRUE))

# Plot F and R qualities
F_read_plots <- plotQualityProfile(fastqFs_path[c(1, 10, 20, 30, 40, 50, 100)])
R_read_plots <- plotQualityProfile(fastqRs_path[c(1, 10, 20, 30, 40, 50, 100)])

# Save to results/microbiome_analysis
# Creat dir if not already present
if (!dir.exists("results/microbiome_analysis")){
  dir.create("results/microbiome_analysis")
}

# Save quality images
ggsave(filename = "results/microbiome_analysis/F_reads_quality.png", 
       plot = F_read_plots)
ggsave(filename = "results/microbiome_analysis/R_reads_quality.png", 
       plot = R_read_plots)


## Filter Reads -------------------------------------------------------------
filter_results <- filterAndTrim(
  fwd = file.path(pathF, fastqFs),
  filt = file.path(filtpathF, fastqFs),
  rev = file.path(pathR, fastqRs),
  filt.rev = file.path(filtpathR, fastqRs),
  truncLen = c(285, 200),
  maxEE = c(2, 2),
  truncQ = 2,
  maxN = 0,
  rm.phix = TRUE,
  compress = TRUE,
  verbose = TRUE,
  multithread = TRUE
)

# Save results to txt file
write.table(x = data.frame(filter_results), 
          file = "results/microbiome_analysis/filter_results.txt", 
          col.names = TRUE)


## Learn the Error Rates -----------------------------------------------------
# Get list of filtered seqeunces
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

# Assign names to filtered sequences
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# Calculate the error rates
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)


## Sample inference and merge of paired-end reads ----------------------------
# Create an empty list of lenghth equal to number of samples
mergers <- vector("list", length(sample_names))

# Change name of each element of list to reflect sample names
names(mergers) <- sample_names

for(sam in sample_names) {
  cat("Processing:", sam, "\n")
  # Remove replicated sequences from F reads
  derepF <- derepFastq(filtFs[[sam]])
  # Remove all sequencing errors to reveal the member of the
  #  sequenced community
  ddF <- dada(derepF, err = errF, multithread = TRUE)
  # Same thing for reverse reads
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err = errR, multithread = TRUE)
  # Merge sequnces
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  # Add merged results to mergers list
  mergers[[sam]] <- merger
}

# Clean temp files
rm(derepF); rm(derepR)

# Create and save sequence table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "results/microbiome_analysis/seqtab.rds")