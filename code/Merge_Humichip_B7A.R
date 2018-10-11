####################################################
# Name: Merge_Humichip_B7A.R
# Author: Ryan Johnson
# Date Created: 8 October, 2018
# Purpose: Take all Humichip data sets for the B7A
#          CHIM study and combine them into a 
#          single data frame
####################################################

library(tidyverse)
library(readxl)

## Read in Humichip Data --------------------------------------------------

# list vector of humichip files
humichip_files <- c(
  "data/raw/Humichip/HumiChipBx123-LTO.txt",
  "data/raw/Humichip/HumiChipBxRpt-LTO.txt"
)

# Create empty data frame
humichip_data <- tibble()

# Start loop that will read in each humichip separately
for (hchip in humichip_files){
  
  # If tibble is empty (first occurence)
  if(is_empty(humichip_data)){
    
    # Read in hchip
    humichip_data <- suppressWarnings(suppressMessages((read_tsv(hchip, guess_max = 100000)))) %>%
      select(-uniqueID,	-proteinGI,	-accessionNo)
  } else {
    
    # Read in hchip and merge into humichip_data
    humichip_temp <-  suppressWarnings(suppressMessages((read_tsv(hchip, guess_max = 100000)))) %>%
      select(-uniqueID,	-proteinGI,	-accessionNo)
    
    humichip_data <- suppressWarnings(suppressMessages((full_join(humichip_data, humichip_temp, 
                               by = c("Genbank ID", "gene", "species", "lineage", "annotation", 
                                      "geneCategory", "subcategory1", "subcategory2")))))
  }
}

# Clean
rm(humichip_temp, hchip, humichip_files)

# Remove any extraneous text at end of sample headers
colnames(humichip_data) <- str_replace(string = colnames(humichip_data), 
                                       pattern = "^[Xx]", 
                                       replacement = "B7A")

colnames(humichip_data) <- str_replace(string = colnames(humichip_data), 
                                       pattern = "(B7A.*)--.*", 
                                       replacement = "\\1")


# Remove duplicates (remove one with least number of probes)
dups <- str_subset(pattern = "^B.*2$", string = colnames(humichip_data)) 
origs <- str_replace(string = dups, pattern = "2$", replacement = "")

humichip_dups_to_remove <- humichip_data %>% # Pull out the dups with the least number positive probes
  select(c(dups, origs)) %>%
  gather(study_ID, value) %>%
  filter(!is.na(value)) %>%
  mutate(dup_group = ifelse(grepl(pattern = "2$", x = study_ID), 2, 1)) %>%
  mutate(study_ID_group = str_replace(string = study_ID, pattern = "2$", replacement = "")) %>%
  group_by(study_ID, dup_group, study_ID_group) %>%
  count() %>%
  ungroup() %>%
  group_by(study_ID_group) %>%
  filter(n == min(n)) %>%
  pull(study_ID)

# Filter the Humichip Data
humichip_data <- humichip_data %>%
  select(-humichip_dups_to_remove)

# Clean
rm(humichip_dups, dups, humichip_dups_to_remove, origs)
