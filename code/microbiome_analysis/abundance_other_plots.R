# abundance_other_plots.R
# Ryan Johnson
# 5 Nov 2019
# Analyze polished microbiome data to produce 
#  relative abundance plots

## Load Packages -------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(ggthemes)
library(RColorBrewer)



## Load Data and set Variables -----------------------------------------------
# Data load
ps2 <- readRDS("results/microbiome_analysis/ps2.rds")

# Variables
num_top_genus <- 10 # How many most abundant genera to consider in plots



## Covert to Relative Abundance ----------------------------------------------
# Combine data to Phylum Level
ps_phylum <- tax_glom(ps2, "Phylum", NArm = TRUE)

# Convert phylum data to Relative Abundacne
ps_phylum_relabun <- transform_sample_counts(ps_phylum, 
                                             function(OTU) OTU/sum(OTU) * 100)
# Convert phylum relative abundance to data frame
taxa_abundance_table_phylum <- psmelt(ps_phylum_relabun)

# Combine data to Genus Level
ps_genus <- tax_glom(ps2, "Genus", NArm = TRUE)

# Convert genus data to relative abundance
ps_genus_relabun <- transform_sample_counts(ps_genus, 
                                            function(OTU) OTU/sum(OTU) * 100)

# Get top N genera
topN_genus <- names(sort(taxa_sums(ps_genus), decreasing=TRUE))[1:num_top_genus]

# Extract the top N genera
ps_genus_topN <- prune_taxa(topN_genus, ps_genus_relabun)

# Convert top genus relative abundance to data frame
taxa_abundance_table_genus <- psmelt(ps_genus_topN)


## Plots --------------------------------------------------------
# Phylum level line plots (mod sev disease)
phylum_line_mod_sev <- taxa_abundance_table_phylum %>% 
  filter(!is.na(Mod.Sev.Diarrhea)) %>% 
  filter(!Study_Day %in% c(28, 36)) %>%
  group_by(Phylum, Study_Day, Mod.Sev.Diarrhea) %>% 
  summarise(mean_abundance = mean(Abundance), 
            sd_abundance = sd(Abundance),
            sem_abundance = sd_abundance / sqrt(n())) %>% 
  ggplot(aes(x = factor(Study_Day), y = mean_abundance, group = Phylum, color = Phylum)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = mean_abundance - sem_abundance,
                    ymax = mean_abundance + sem_abundance),
                width = 0.2) +
  labs(title = "Phylum Abundance", 
       subtitle = "Y: Mod-Sev Disease\nN: Mild Disease",
       caption = "Error bars: SEM",
       y = "Mean Relative Abundance", 
       x = "Study Day") +
  facet_grid(~Mod.Sev.Diarrhea) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(plot = phylum_line_mod_sev, 
       filename = "results/microbiome_analysis/abundance_other_plots/phylum_line_mod_sev.png",
       height = 6, 
       width = 7)

########################################################

# Phylum level line plots (naive)
phylum_line_naive <- taxa_abundance_table_phylum %>% 
  filter(!is.na(Naïve)) %>% 
  filter(!Study_Day %in% c(28, 36)) %>%
  group_by(Phylum, Study_Day, Naïve) %>% 
  summarise(mean_abundance = mean(Abundance), 
            sd_abundance = sd(Abundance),
            sem_abundance = sd_abundance / sqrt(n())) %>% 
  ggplot(aes(x = factor(Study_Day), y = mean_abundance, group = Phylum, color = Phylum)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = mean_abundance - sem_abundance,
                    ymax = mean_abundance + sem_abundance),
                width = 0.2) +
  labs(title = "Phylum Abundance", 
       subtitle = "Y: Naive\nN: Not-Naive",
       caption = "Error bars: SEM",
       y = "Mean Relative Abundance", 
       x = "Study Day") +
  facet_grid(~Naïve) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(plot = phylum_line_naive, 
       filename = "results/microbiome_analysis/abundance_other_plots/phylum_line_naive.png",
       height = 6, 
       width = 7)

#######################################################

# Phylum level line plots (abx)
taxa_abundance_table_phylum$pre_post_abx <- factor(taxa_abundance_table_phylum$pre_post_abx,
                                                   levels = c("Pre", "Post"))

phylum_line_abx <- taxa_abundance_table_phylum %>% 
  filter(!is.na(pre_post_abx)) %>% 
  filter(!Study_Day %in% c(28, 36)) %>%
  group_by(Phylum, Study_Day, pre_post_abx) %>% 
  summarise(mean_abundance = mean(Abundance), 
            sd_abundance = sd(Abundance),
            sem_abundance = sd_abundance / sqrt(n())) %>% 
  ggplot(aes(x = factor(Study_Day), y = mean_abundance, group = Phylum, color = Phylum)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = mean_abundance - sem_abundance,
                    ymax = mean_abundance + sem_abundance),
                width = 0.2) +
  labs(title = "Phylum Abundance", 
       subtitle = "Pre or Post Abx Treatment",
       caption = "Error bars: SEM",
       y = "Mean Relative Abundance", 
       x = "Study Day") +
  facet_grid(~pre_post_abx) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(plot = phylum_line_abx, 
       filename = "results/microbiome_analysis/abundance_other_plots/phylum_line_abx.png",
       height = 6, 
       width = 7)

####################################################

# genus level line plots (mod sev disease)
genus_line_mod_sev <- taxa_abundance_table_genus %>% 
  filter(!is.na(Mod.Sev.Diarrhea)) %>% 
  filter(!Study_Day %in% c(28, 36)) %>%
  group_by(Genus, Study_Day, Mod.Sev.Diarrhea) %>% 
  summarise(mean_abundance = mean(Abundance), 
            sd_abundance = sd(Abundance),
            sem_abundance = sd_abundance / sqrt(n())) %>% 
  ggplot(aes(x = factor(Study_Day), y = mean_abundance, group = Genus, color = Genus)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = mean_abundance - sem_abundance,
                    ymax = mean_abundance + sem_abundance),
                width = 0.2) +
  labs(title = paste0("Genus Top ", num_top_genus, " Relative Abundance"), 
       subtitle = "Y: Mod-Sev Disease\nN: Mild Disease",
       caption = "Error bars: SEM",
       y = "Mean Relative Abundance", 
       x = "Study Day") +
  facet_grid(~Mod.Sev.Diarrhea) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(plot = genus_line_mod_sev, 
       filename = "results/microbiome_analysis/abundance_other_plots/genus_line_mod_sev.png",
       height = 6, 
       width = 7)

########################################################

# Genus level line plots (naive)
genus_line_naive <- taxa_abundance_table_genus %>% 
  filter(!is.na(Naïve)) %>% 
  filter(!Study_Day %in% c(28, 36)) %>%
  group_by(Genus, Study_Day, Naïve) %>% 
  summarise(mean_abundance = mean(Abundance), 
            sd_abundance = sd(Abundance),
            sem_abundance = sd_abundance / sqrt(n())) %>% 
  ggplot(aes(x = factor(Study_Day), y = mean_abundance, group = Genus, color = Genus)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = mean_abundance - sem_abundance,
                    ymax = mean_abundance + sem_abundance),
                width = 0.2) +
  labs(title = paste0("Genus Top ", num_top_genus, " Relative Abundance"), 
       subtitle = "Y: Naive\nN: Not-Naive",
       caption = "Error bars: SEM",
       y = "Mean Relative Abundance", 
       x = "Study Day") +
  facet_grid(~Naïve) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(plot = genus_line_naive, 
       filename = "results/microbiome_analysis/abundance_other_plots/genus_line_naive.png",
       height = 6, 
       width = 7)

#######################################################

# Genus level line plots (abx)
taxa_abundance_table_genus$pre_post_abx <- factor(taxa_abundance_table_genus$pre_post_abx,
                                                   levels = c("Pre", "Post"))

genus_line_abx <- taxa_abundance_table_genus %>% 
  filter(!is.na(pre_post_abx)) %>% 
  filter(!Study_Day %in% c(28, 36)) %>%
  group_by(Genus, Study_Day, pre_post_abx) %>% 
  summarise(mean_abundance = mean(Abundance), 
            sd_abundance = sd(Abundance),
            sem_abundance = sd_abundance / sqrt(n())) %>% 
  ggplot(aes(x = factor(Study_Day), y = mean_abundance, group = Genus, color = Genus)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = mean_abundance - sem_abundance,
                    ymax = mean_abundance + sem_abundance),
                width = 0.2) +
  labs(title = paste0("Genus Top ", num_top_genus, " Relative Abundance"), 
       subtitle = "Pre or Post Abx Treatment",
       caption = "Error bars: SEM",
       y = "Mean Relative Abundance", 
       x = "Study Day") +
  facet_grid(~pre_post_abx) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(plot = genus_line_abx, 
       filename = "results/microbiome_analysis/abundance_other_plots/genus_line_abx.png",
       height = 6, 
       width = 7)
