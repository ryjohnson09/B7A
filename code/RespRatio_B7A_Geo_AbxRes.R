################################################
# Name: RespRatio_B7A_Geo_AbxRes.R
# Author: Ryan Johnson
# Date Created: 8 February, 2019
# Purpose: Generate RR from geochip B7A data for
#  just the antibiotic/drug resistance probes
################################################

library(tidyverse)

## Read in Data -------------------------------------------------------------------
geochip <- suppressWarnings(suppressMessages(read_tsv("data/processed/Merged_geochip_B7A.tsv")))
metadata <- suppressWarnings(suppressMessages(read_tsv("data/processed/B7A_metadata.tsv")))

## Filter metadata for samples in geochip ----------------------------------------
metadata <- metadata %>%
  filter(glomics_ID %in% colnames(geochip))

## Only keep matched isolates -----------------------------------------------------
# I am most interested in comparing the B samples (acute illness right before abx)
# to C samples
metadata_matched <- metadata %>% 
  filter(Sample %in% c("B", "C")) %>% 
  group_by(`Subject ID`) %>% 
  filter(n() == 2) %>% 
  ungroup()

rm(metadata)

## Filter the geochip data to include only samples in metadata_matched --------------------------------
geochip_filtered <- geochip %>% 
  select_if(colnames(.) %in% c("Genbank.ID", "Gene", "Organism", "Lineage",
                               "Gene_category", "Subcategory1",
                               "Subcategory2", metadata_matched$glomics_ID))

rm(geochip)


## Convert the geochip data to relative abundance -------------------------------------------------
geo_relabun <- geochip_filtered %>%
  # only consider probes with a gene category designation
  filter(!is.na(Gene_category)) %>%
  # Convert all values to 1 and 0
  mutate_at(vars(starts_with("B7A")), funs(ifelse(is.na(.), 0, 1))) %>% 
  # Convert all values to relative abundance
  mutate_at(vars(starts_with("B7A")), funs((./sum(.)) * 100))

rm(geochip_filtered)


## Select probes associated with antibioic resistance/drug resistance
geo_abx <- geo_relabun %>% 
  filter(Subcategory1 %in% c("ANTIBIOTIC_RESISTANCE", "DRUG_RESISTANCE"))

rm(geo_relabun)

## Summarize Data ----------------------------------------------------------------------------------
geo_grouped <- geo_abx %>%
  # Select grouping column of interest and remove rest
  select(Gene, starts_with("B7A")) %>%
  # Make long
  gather(key = glomics_ID, value = rel_abun_value, -Gene) %>%
  # Group by category of interest
  group_by(glomics_ID, Gene) %>%
  # Calculate total relative abundance for each category
  summarise(category_abundance = sum(rel_abun_value)) %>%
  # Add in metadata
  left_join(., metadata_matched, by = "glomics_ID") %>% 
  ungroup()

rm(geo_relabun)

## Calculate Response Ratio ------------------------------------------------------------------------

geo_RR <- geo_grouped %>% 
  # Group by category and visit
  group_by(Gene, Sample) %>% 
  summarise(mean_signal = mean(category_abundance),
            sd_signal = sd(category_abundance),
            n = sum(!is.na(category_abundance))) %>% 
  
  # Spread the signal mean by visit number
  ungroup() %>%
  group_by(Gene) %>%
  spread(Sample, mean_signal) %>% 
  
  # Rename mean columns
  rename(group1_mean = colnames(.)[length(colnames(.)) - 1],
         group2_mean = colnames(.)[length(colnames(.))]) %>% 
  
  # Spread the sd and n columns by visit
  mutate(sd_group1 = ifelse(!is.na(group1_mean), sd_signal, NA)) %>%
  mutate(sd_group2 = ifelse(!is.na(group2_mean), sd_signal, NA)) %>%
  mutate(n_group1 = ifelse(!is.na(group1_mean), n, NA)) %>%
  mutate(n_group2 = ifelse(!is.na(group2_mean), n, NA)) %>%
  select(-sd_signal, -n) %>% 
  
  # Compress NAs
  ungroup() %>%
  group_by(Gene) %>%
  summarise_all(funs(sum(., na.rm = T))) %>% 
  
  # Must have at least __ observations in each subcategory
  filter(n_group1 >= 10) %>%
  filter(n_group2 >= 10) %>%
  
  # Calculate SEM for each mean
  mutate(SEM_group1 = sd_group1 / sqrt(n_group1)) %>%
  mutate(SEM_group2 = sd_group2 / sqrt(n_group2)) %>%
  
  # Calculate the Response Ratio (RR)
  mutate(RR = log(group2_mean / group1_mean)) %>%
  
  # Calculate the Standard error for the RR
  mutate(SE_RR = sqrt((SEM_group1**2 / group1_mean**2) + (SEM_group2**2 / group2_mean**2))) %>%
  
  # Calcualte the 95% confidence interval for each RR
  mutate(CI95 = abs(1.96 * SE_RR)) %>%
  
  # Add in keeper column if does not overlap 0
  mutate(keeper = ifelse(0 > (RR - CI95) & 0 < (RR + CI95), "No", "Yes")) %>%
  
  # Make labels pretty
  mutate(pretty_cat = str_to_title(Gene)) %>%
  mutate(pretty_cat = str_replace_all(pretty_cat, "_"," ")) %>%
  
  # Factor columns
  mutate(pretty_cat = fct_reorder(pretty_cat, RR)) %>%
  ungroup()

## Remove any cateogories where the mean was below threshold ----------------
geo_RR_filtered <- geo_RR %>% 
  filter(group1_mean > 0.01) %>% 
  filter(group2_mean > 0.01)


## Plot ---------------------------------------------------------------------

RR_plot <- ggplot(data = geo_RR_filtered) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  
  # points and error bar
  geom_point(aes(x = pretty_cat, y = RR), 
             size = 4) +
  geom_errorbar(aes(ymin = RR - CI95, 
                    ymax = RR + CI95, 
                    x = pretty_cat),
                width = 0.25) +
  
  # Group labels
  annotate(geom = "text", label = "Sample B", x = Inf, y = -Inf, hjust = 0, vjust = 1, 
           size = 5, color = "red", fontface = 2) +
  annotate(geom = "text", label = "Sample C", x = Inf, y = Inf, hjust = 1, vjust = 1, 
           size = 5, color = "red", fontface = 2) +
  
  # plot labels
  labs(title = "Response Ratio",
       x = "Category",
       y = "Response Ratio") +
  
  theme_minimal() +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5)
  )

RR_plot

ggsave(plot = RR_plot, filename = "results/figures/RespRatio_Geo_AbxRes_gene_B7A.png", height = 8, width = 7)
