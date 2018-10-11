####################################################
# Name: Merge_Geochip_B7A.R
# Author: Ryan Johnson
# Date Created: 8 October, 2018
# Purpose: Take all Geochip data sets for the B7A
#          CHIM study and combine them into a 
#          single data frame
####################################################

library(tidyverse)
library(readxl)

## Read in Geochip Data --------------------------------------------------

# list vector of geochip files
geochip_files <- c(
  "data/raw/Geochip/GeoChipBx1-LTO.txt",
  "data/raw/Geochip/GeoChipBx2-LTO.txt",
  "data/raw/Geochip/GeoChipBx3-LTO.txt",
  "data/raw/Geochip/GeoChipBxRpt-LTO.txt"
)

# Create empty data frame
geochip_data <- tibble()

# Start loop that will read in each geochip separately
for (gchip in geochip_files){
  
  # If tibble is empty (first occurence)
  if(is_empty(geochip_data)){
    
    # Read in gchip
    geochip_data <- suppressWarnings(suppressMessages((read_tsv(gchip, guess_max = 100000))))

  } else {
    
    # Read in gchip and merge into geochip_data
    geochip_temp <-  suppressWarnings(suppressMessages((read_tsv(gchip, guess_max = 100000))))
    
    geochip_data <- suppressWarnings(suppressMessages((full_join(geochip_data, geochip_temp, 
                                                                  by = c("Genbank ID", "Gene", "Organism", "Gene_category", 
                                                                         "Subcategory1", "Subcategory2", "Lineage")))))
  }
}

# Clean
rm(geochip_temp, gchip, geochip_files)


# Remove any extraneous text at end of sample headers
colnames(geochip_data) <- str_replace(string = colnames(geochip_data), 
                                       pattern = "^[Xx]", 
                                       replacement = "B7A")

colnames(geochip_data) <- str_replace(string = colnames(geochip_data), 
                                       pattern = "(B7A.*)--.*", 
                                       replacement = "\\1")

colnames(geochip_data) <- str_replace(string = colnames(geochip_data), # remove internal "x"
                                      pattern = "(B7A.*)x(.*)",
                                      replacement = "\\1\\2")



# Remove duplicates (remove one with least number of probes)
dups <- str_subset(pattern = "^B.*2$", string = colnames(geochip_data)) 
origs <- str_replace(string = dups, pattern = "2$", replacement = "")

geochip_dups_to_remove <- geochip_data %>% # Pull out the dups with the least number positive probes
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

# Filter the geochip Data
geochip_data <- geochip_data %>%
  select(-geochip_dups_to_remove)

# Remove the "2" from any remaining sample names
colnames(geochip_data) <- str_replace(string = colnames(geochip_data), 
                                       pattern = "(B7A.*)2$", 
                                       replacement = "\\1")

# Clean
rm(dups, origs, geochip_dups_to_remove)


# Fix Gene_Categories
geochip_data$Gene_category <- toupper(geochip_data$Gene_category)
geochip_data$Gene_category <- gsub(pattern = " ", replacement = "_", geochip_data$Gene_category)


# Fix Subcategory1
geochip_data$Subcategory1 <- toupper(geochip_data$Subcategory1)
geochip_data$Subcategory1 <- gsub(pattern = " ", replacement = "_", geochip_data$Subcategory1)

geochip_data <- geochip_data %>%
  mutate(Subcategory1 = ifelse(Subcategory1 == "EFFECTOR", "EFFECTOR_PROTEIN", Subcategory1))

# Fix Subcategory2
geochip_data$Subcategory2 <- toupper(geochip_data$Subcategory2)
geochip_data$Subcategory2 <- gsub(pattern = " ", replacement = "_", geochip_data$Subcategory2)

geochip_data <- geochip_data %>%
  mutate(Subcategory2 = ifelse(Subcategory2 == "OOMYCETES", "OOMYCETE", Subcategory2)) %>%
  mutate(Subcategory2 = ifelse(Subcategory2 == "TRANSPORTER", "TRANSPORT", Subcategory2))


## Fix Sample names (must be B7A###[ABC] format)
geochip_data <- geochip_data %>%
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
write_tsv(geochip_data, "data/processed/Merged_geochip_B7A.tsv")
