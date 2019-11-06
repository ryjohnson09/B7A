# alpha_diversity_analysis.R
# Ryan Johnson
# 6 Nov 2019
# Alpha diversity metrics/plots for B7A
#  microbiome data

## Load Libraries ---------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(ggthemes)
theme_set(theme_bw())


## Load Data --------------------------------------------------------------------
ps2 <- readRDS("results/microbiome_analysis/ps2.rds")
metadata <- read_csv("data/processed/microbiome_analysis/B7A_metadata_clean.csv")


## Calculate Alpha Diversity Metrics --------------------------------------------
alpha_diversity <- estimate_richness(ps2, 
                                     split = TRUE, 
                                     measures = c("Observed", "Shannon", "InvSimpson")) %>% 
  as_tibble(rownames = "TUBE") %>% 
  mutate(TUBE = str_replace(TUBE, "\\.S", "_S")) %>% 
  mutate(TUBE = str_replace(TUBE, "\\.", "-"))

# Merge with metadata
alpha_diversity_metadata <- alpha_diversity %>% 
  left_join(metadata, by = c("TUBE"))


## Plots ---[observed ASVs]------------------------------------------------------
# Observed ASVs over time, mod sev diarrhea
alpha_diversity_metadata$`Mod-Sev Diarrhea` <- 
  factor(alpha_diversity_metadata$`Mod-Sev Diarrhea`, levels = c("Y", "N"))

observed_mod_sev <- alpha_diversity_metadata %>% 
  ggplot(aes(x = factor(Study_Day), y = Observed, fill = `Mod-Sev Diarrhea`)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
              pch = 21, color = "black") +
  labs(x = "Study Day", 
       y = "Observed ASVs",
       title = "Observed ASVs",
       subtitle = "Mod-Severe Disease (Yes or No)",
       fill = "Mod to Severe\nDiarrhea?") +
  scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(plot = observed_mod_sev,
       filename = "results/microbiome_analysis/alpha_diversity_analysis/observed_mod_sev.png",
       width = 7,
       height = 5)

#################################################

# Observed ASVs over time, naive
alpha_diversity_metadata$Naïve <- factor(alpha_diversity_metadata$Naïve, levels = c("Y", "N"))

observed_naive <- alpha_diversity_metadata %>% 
  ggplot(aes(x = factor(Study_Day), y = Observed, fill = Naïve)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
              pch = 21, color = "black") +
  labs(x = "Study Day", 
       y = "Observed ASVs",
       title = "Observed ASVs",
       subtitle = "Naïve? (Yes or No)",
       fill = "Naïve?") +
  scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(plot = observed_naive,
       filename = "results/microbiome_analysis/alpha_diversity_analysis/observed_naive.png",
       width = 7,
       height = 5)

#################################################

# Observed ASVs over time, abx
alpha_diversity_metadata$pre_post_abx <- 
  factor(alpha_diversity_metadata$pre_post_abx, levels = c("Pre", "Post"))

observed_abx <- alpha_diversity_metadata %>% 
  ggplot(aes(x = factor(Study_Day), y = Observed, fill = pre_post_abx)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
              pch = 21, color = "black") +
  labs(x = "Study Day", 
       y = "Observed ASVs",
       title = "Observed ASVs",
       subtitle = "Pre or Post Abx Treatment",
       fill = "Pre or Post\nAntibiotics") +
  scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(plot = observed_abx,
       filename = "results/microbiome_analysis/alpha_diversity_analysis/observed_abx.png",
       width = 7,
       height = 5)



## Plots ---[InvSimpson]---------------------------------------------------
# InvSimpson over time, mod sev
invSimpson_mod_sev <- alpha_diversity_metadata %>% 
  ggplot(aes(x = factor(Study_Day), y = InvSimpson, fill = `Mod-Sev Diarrhea`)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
              pch = 21, color = "black") +
  labs(x = "Study Day", 
       y = "InvSimpson",
       title = "InvSimpson",
       subtitle = "Mod-Severe Disease (Yes or No)",
       fill = "Mod to Severe\nDiarrhea?") +
  scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(plot = invSimpson_mod_sev,
       filename = "results/microbiome_analysis/alpha_diversity_analysis/invSimpson_mod_sev.png",
       width = 7,
       height = 5)

#################################################

# InvSimpson over time, naive
invSimpson_naive <- alpha_diversity_metadata %>% 
  ggplot(aes(x = factor(Study_Day), y = InvSimpson, fill = Naïve)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
              pch = 21, color = "black") +
  labs(x = "Study Day", 
       y = "InvSimpson",
       title = "InvSimpson",
       subtitle = "Naïve? (Yes or No)",
       fill = "Naïve?") +
  scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(plot = invSimpson_naive,
       filename = "results/microbiome_analysis/alpha_diversity_analysis/invSimpson_naive.png",
       width = 7,
       height = 5)

#################################################

# InvSimpson over time, abx
invSimpson_abx <- alpha_diversity_metadata %>% 
  ggplot(aes(x = factor(Study_Day), y = InvSimpson, fill = pre_post_abx)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
              pch = 21, color = "black") +
  labs(x = "Study Day", 
       y = "InvSimpson",
       title = "InvSimpson",
       subtitle = "Pre or Post Abx Treatment",
       fill = "Pre or Post\nAntibiotics") +
  scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(plot = invSimpson_abx,
       filename = "results/microbiome_analysis/alpha_diversity_analysis/invSimpson_abx.png",
       width = 7,
       height = 5)

## Plots ---[Shannon]---------------------------------------------------
# Shannon over time, mod sev
shannon_mod_sev <- alpha_diversity_metadata %>% 
  ggplot(aes(x = factor(Study_Day), y = Shannon, fill = `Mod-Sev Diarrhea`)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
              pch = 21, color = "black") +
  labs(x = "Study Day", 
       y = "Shannon",
       title = "Shannon",
       subtitle = "Mod-Severe Disease (Yes or No)",
       fill = "Mod to Severe\nDiarrhea?") +
  scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(plot = shannon_mod_sev,
       filename = "results/microbiome_analysis/alpha_diversity_analysis/shannon_mod_sev.png",
       width = 7,
       height = 5)

#################################################

# Shannon over time, naive
shannon_naive <- alpha_diversity_metadata %>% 
  ggplot(aes(x = factor(Study_Day), y = Shannon, fill = Naïve)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
              pch = 21, color = "black") +
  labs(x = "Study Day", 
       y = "Shannon",
       title = "Shannon",
       subtitle = "Naïve? (Yes or No)",
       fill = "Naïve?") +
  scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(plot = shannon_naive,
       filename = "results/microbiome_analysis/alpha_diversity_analysis/shannon_naive.png",
       width = 7,
       height = 5)

#################################################

# Shannon over time, abx
shannon_abx <- alpha_diversity_metadata %>% 
  ggplot(aes(x = factor(Study_Day), y = Shannon, fill = pre_post_abx)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
              pch = 21, color = "black") +
  labs(x = "Study Day", 
       y = "Shannon",
       title = "Shannon",
       subtitle = "Pre or Post Abx Treatment",
       fill = "Pre or Post\nAntibiotics") +
  scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(plot = shannon_abx,
       filename = "results/microbiome_analysis/alpha_diversity_analysis/shannon_abx.png",
       width = 7,
       height = 5)
