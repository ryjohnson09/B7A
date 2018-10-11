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

## Fix times ---------------------------------------------------

B7A_fixed_times <- B7A_metadata %>%
  
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

