## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(rstatix)
library(ggpubr)
library(devtools)

## source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

## load data
dat <- read_csv("Data/Tissue_paper/cSCC_miRNA_irradiation.csv") %>% 
  mutate(treatment = factor(treatment, levels = c("unirradiated", "irradiated"))) %>%
  mutate(cluster = factor(cluster, levels = c(1,4,3), labels = c("Cluster 1", "Cluster 4", "Additional miRNAs")))

## stat test
dat %>% group_by(cluster) %>%
  anova_test(log2(expression)~treatment*miRNA) 

## plot data
svg("Results/Tissue_paper/cSCC_miRNA_irradiation.svg", width=9, height=6)
dat %>% 
  ggplot(aes(miRNA, log2(expression), fill = treatment)) + 
  geom_boxplot(width = 0.6, outlier.shape = NA, position = position_dodge())+
  geom_jitter(alpha = 0.5, position = position_dodge(0.6))+
  facet_grid (.~ cluster, scales = "free", space = "free_x") +
  scale_fill_manual(values = c("white", "grey40")) +
  theme_bw() +
  # theme_PhD(axis.text.size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  scale_y_continuous(guide = "axis_minor",limits = c(3,12.53), breaks = seq(4,12,2)) +
  ylab("log2 miRNA expression")
dev.off()



