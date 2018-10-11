################################################
# Name: Process_B7A_Metadata.R
# Author: Ryan Johnson
# Date Created: 11 October, 2018
# Purpose: Manipulate the B7A metadata so that
#   it works nicely with Humichip and Geochip 
#   data
################################################

library(tidyverse)
library(readxl)
library(lubridate)

## Read in data ------------------------------------------------
B7A_metadata <- read_excel("data/raw/B7A_Listing_(Glomics)_10OCT2018.xlsx")

humichip_samples <- read_tsv("data/processed/Merged_humichip_B7A.tsv", n_max = 1) %>%
  select(starts_with("B7A")) %>%
  colnames(.)

geochip_samples <- read_tsv("data/processed/Merged_geochip_B7A.tsv", n_max = 1) %>%
  select(starts_with("B7A")) %>%
  colnames(.)

## Fix times ---------------------------------------------------

B7A_metadata <- B7A_metadata %>%
  
  # Sample date/time
  mutate(`Sample Time` = format(`Sample Time`, "%H:%M")) %>%
  mutate(`Sample Date` = ymd_hm(paste(`Sample Date`, `Sample Time`))) %>%
  select(-`Sample Time`) %>%
  
  # Inoculate date/time
  mutate(`Time inoculated` = format(`Time inoculated`, "%H:%M")) %>%
  mutate(`Date inoculated` = ymd_hm(paste(`Date inoculated`, `Time inoculated`))) %>%
  select(-`Time inoculated`) %>%
  
  # Antibiotics date/time
  mutate(`Time antibiotics initiated` = format(`Time antibiotics initiated`, "%H:%M")) %>%
  mutate(`Date antibiotics initiated` = ymd_hm(paste(`Date antibiotics initiated`, `Time antibiotics initiated`))) %>%
  select(-`Time antibiotics initiated`)


## Fix Inoculation doses ----------------------------------------
B7A_metadata <- B7A_metadata %>%
  
  # Target dose
  mutate(`Inoculum dose (Target)` = as.numeric(str_replace(
    string = `Inoculum dose (Target)`, 
    pattern = "\\sx\\s10", 
    replacement = "e+"))) %>%
  
  # Target dose
  mutate(`Inoculum dose (Actual)` = as.numeric(str_replace(
    string = `Inoculum dose (Actual)`, 
    pattern = "\\sx\\s10", 
    replacement = "e+")))


# Ensure all humichip samples & Geochip have a matching ID in metadata
ID_test <- 1
if(!is_empty(setdiff(humichip_samples, paste0(B7A_metadata$`Subject ID`, B7A_metadata$Sample)))){
  ID_test <- 2
}
if(!is_empty(setdiff(geochip_samples, paste0(B7A_metadata$`Subject ID`, B7A_metadata$Sample)))){
  ID_test <- 3
}

# Write final version
if(ID_test == 1){
  write_tsv(x = B7A_metadata, path = "data/processed/B7A_metadata.tsv")
} else if (ID_test == 2){
  message("Humichip sample not in B7A metadata excel file")
} else {
  message("Geochip sample not in B7A metadata excel file")
}

