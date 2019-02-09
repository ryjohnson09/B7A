################################################
# Name: RespRatio_B7A_Humi.R
# Author: Ryan Johnson
# Date Created: 8 February, 2019
# Purpose: Generate RR from humichip B7A data
################################################

library(tidyverse)

## Read in Data -------------------------------------------------------------------
humichip <- suppressWarnings(suppressMessages(read_tsv("data/processed/Merged_humichip_B7A.tsv")))
metadata <- suppressWarnings(suppressMessages(read_tsv("data/processed/B7A_metadata.tsv")))

## Filter metadata for samples in Humichip ----------------------------------------
metadata <- metadata %>%
  filter(glomics_ID %in% colnames(humichip))

## Only keep matched isolates -----------------------------------------------------
# I am most interested in comparing the B samples (acute illness right before abx)
# to C samples
metadata_matched <- metadata %>% 
  filter(Sample %in% c("B", "C")) %>% 
  group_by(`Subject ID`) %>% 
  filter(n() == 2) %>% 
  ungroup()

rm(metadata)

## Filter the humichip data to include only samples in metadata_matched --------------------------------
humichip_filtered <- humichip %>% 
  select_if(colnames(.) %in% c("Genbank.ID", "gene", "species", "lineage",
                               "annotation", "geneCategory", "subcategory1",
                               "subcategory2", metadata_matched$glomics_ID))

rm(humichip)


## Convert the Humichip data to relative abundance -------------------------------------------------
humi_relabun <- humichip_filtered %>%
  # only consider probes with a gene category designation
  filter(!is.na(geneCategory)) %>%
  # Convert all values to 1 and 0
  mutate_at(vars(starts_with("B7A")), funs(ifelse(is.na(.), 0, 1))) %>% 
  # Convert all values to relative abundance
  mutate_at(vars(starts_with("B7A")), funs((./sum(.)) * 100))

rm(humichip_filtered)


## Summarize Data ----------------------------------------------------------------------------------
humi_grouped <- humi_relabun %>%
  # Select grouping column of interest and remove rest
  select(geneCategory, starts_with("B7A")) %>%
  # Make long
  gather(key = glomics_ID, value = rel_abun_value, -geneCategory) %>%
  # Group by category of interest
  group_by(glomics_ID, geneCategory) %>%
  # Calculate total relative abundance for each category
  summarise(category_abundance = sum(rel_abun_value)) %>%
  # Add in metadata
  left_join(., metadata_matched, by = "glomics_ID") %>% 
  ungroup()

rm(humi_relabun)

## Calculate Response Ratio ------------------------------------------------------------------------

humi_RR <- humi_grouped %>% 
  # Group by category and visit
  group_by(geneCategory, Sample) %>% 
  summarise(mean_signal = mean(category_abundance),
            sd_signal = sd(category_abundance),
            n = sum(!is.na(category_abundance))) %>% 
  
  # Spread the signal mean by visit number
  ungroup() %>%
  group_by(geneCategory) %>%
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
  group_by(geneCategory) %>%
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
  mutate(pretty_cat = str_to_title(geneCategory)) %>%
  mutate(pretty_cat = str_replace_all(pretty_cat, "_"," ")) %>%
  
  # Factor columns
  mutate(pretty_cat = fct_reorder(pretty_cat, RR)) %>%
  ungroup()

## Remove any cateogories where the mean was below threshold ----------------
humi_RR_filtered <- humi_RR %>% 
  filter(group1_mean > 0.01) %>% 
  filter(group2_mean > 0.01)


## Plot ---------------------------------------------------------------------

RR_plot <- ggplot(data = humi_RR_filtered) +
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

ggsave(plot = RR_plot, filename = "results/figures/RespRatio_Humi_geneCat_B7A.png", height = 8, width = 7)
