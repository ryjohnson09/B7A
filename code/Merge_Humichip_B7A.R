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

# Remove the "2" from any remaining sample names
colnames(humichip_data) <- str_replace(string = colnames(humichip_data), 
                                       pattern = "(B7A.*)2$", 
                                       replacement = "\\1")


# Clean
rm(dups, origs, humichip_dups_to_remove)



## Fixing Category names ------------------------------------------------------

# geneCategory
humichip_data$geneCategory <- toupper(humichip_data$geneCategory)
humichip_data$geneCategory <- gsub(pattern = " ", replacement = "_", humichip_data$geneCategory)

# annotation
humichip_data$annotation <- toupper(humichip_data$annotation)
humichip_data$annotation <- gsub(pattern = " ", replacement = "_", humichip_data$annotation)

# subcategory1
humichip_data$subcategory1 <- toupper(humichip_data$subcategory1)
humichip_data$subcategory1 <- gsub(pattern = " ", replacement = "_", humichip_data$subcategory1)

# subcategory2
humichip_data$subcategory2 <- toupper(humichip_data$subcategory2)
humichip_data$subcategory2 <- gsub(pattern = " ", replacement = "_", humichip_data$subcategory2)



## Fix Sample names (must be B7A###[ABC] format)
humichip_data <- humichip_data %>%
  rename(B7A037A = B7A37A,
         B7A053A = B7A53A,
         B7A065A = B7A65A,
         B7A066A = B7A66A,
         B7A004B = B7A4B,
         B7A050B = B7A50B,
         B7A070B = B7A70B,
         B7A006C = B7A6C,
         B7A010C = B7A10C,
         B7A012C = B7A12C,
         B7A043C = B7A43C)
  
  
  
# Return compiled data frame
write_tsv(humichip_data, "data/processed/Merged_humichip_B7A.tsv")
