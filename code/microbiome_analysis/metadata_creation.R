# metadata_reation.R
# Ryan Johnson
# 4 Nov 2019
# Generate Metadata to analyze the B7A micrbiome data

# Load Libraries ----------------------------------------------------------
library(tidyverse)
library(readxl)
library(lubridate)

## Read in Data and Reformat ----------------------------------------------
metadata <- read_xlsx("data/raw/microbiome_sequencing/ED_B7A01_Microbiome_MergedWithClinicalData_2019.01.11.xlsx") %>% 
  select(TUBE, `Study ID`, Sample_Date, Sample_Time, Study_Day, Max_24num, 
         Max_24vol, DiseaseSeverityScore, `Mod-Sev Diarrhea`, Naïve) %>% 
  # Sample date/time
  mutate(`Sample_Time` = format(`Sample_Time`, "%H:%M")) %>%
  mutate(`Sample_Date` = ymd_hm(paste(`Sample_Date`, `Sample_Time`))) %>%
  select(-`Sample_Time`) %>% 
  # Study ID fix
  mutate(Study_ID = paste0("B7A", `Study ID`)) %>% 
  select(-`Study ID`) %>% 
  # Reorder
  select(TUBE, Study_ID, everything())

# Read in additional samples (obtaioned from chad later)
metadata2 <- read_xlsx("data/raw/microbiome_sequencing/MissingSamples_ClinicalData_2019.10.07.xlsx") %>% 
  select(TUBE, `Study ID`, Sample_Date, Sample_Time, Study_Day, Max_24num, 
         Max_24vol, DiseaseSeverityScore, `Mod-Sev Diarrhea`, Naïve) %>% 
  # Sample date/time
  mutate(`Sample_Time` = format(`Sample_Time`, "%H:%M")) %>%
  mutate(`Sample_Date` = mdy_hm(paste(`Sample_Date`, `Sample_Time`))) %>%
  select(-`Sample_Time`) %>% 
  # Study ID fix
  mutate(Study_ID = paste0("B7A", `Study ID`)) %>% 
  select(-`Study ID`) %>% 
  # Reorder
  select(TUBE, Study_ID, everything())

# Merge the two metadata files:
metadata_merge <- metadata %>% 
  full_join(metadata2)


## Add in antibiotic/innoculation date admin time --------------------------------------------
abx_time <- read_xlsx("data/raw/microbiome_sequencing/B7A AbxTimes_09OCT2019.xlsx") %>% 
  mutate(abx_time = format(EAA1DRTM, "%H:%M")) %>% 
  mutate(abx_date = ymd_hm(paste(EAA1DRDT, abx_time))) %>% 
  select(-EAA1DRTM, -EAA1DRDT, -abx_time) %>% 
  filter(COHORT == 2) %>% 
  select(-COHORT)

# Merge abx time to metadata:
metadata_merge_abx <- metadata_merge %>% 
  left_join(abx_time, by = c("Study_ID" = "PATID")) %>% 
  select(TUBE, Study_ID, Sample_Date, Study_Day, abx_date, everything())

# Time since/before abx treatment
metadata_merge_abx <- metadata_merge_abx %>% 
  mutate(pre_post_abx = ifelse(Sample_Date < abx_date, "Pre", "Post")) %>% 
  mutate(hours_from_abx = as.numeric(interval(abx_date, Sample_Date), unit = "hours"))

# Time since innoculation
metadata_merge_abx <- metadata_merge_abx %>% 
  mutate(innoc_time = ymd_hm(c("2016-06-14-10-23"))) %>% 
  mutate(hours_from_innoc = as.numeric(interval(innoc_time, Sample_Date), unit = "hours"))


# Write to file -----------------------------------------------------------------------------
write_csv(metadata_merge_abx, "data/processed/microbiome_analysis/B7A_metadata_clean.csv")
