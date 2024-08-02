## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(rstatix)
library(ggpubr)
library(devtools)

## source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

## load data
dat <- read_csv("Data/Tissue_paper/cSCC_mRNA_irradiation.csv") %>% 
  mutate(treatment = factor(treatment, levels = c("unirradiated", "irradiated")))

## calculate ratio of EMT markers
EMT_markers <- dat %>% 
  filter(gene == "CDH1") %>%
  group_by(treatment) %>%
  summarize(CDH1_mean = mean(dCT)) %>% 
  inner_join(filter(dat, gene == "CDH2")) %>%
  mutate(ratio = dCT - CDH1_mean)
  

## convert dCT to ddCT
ddCT <- dat %>% 
  filter(treatment == "unirradiated") %>%
  group_by(gene) %>%
  summarize(ctrl_dCT = mean(dCT, na.rm = TRUE)) %>%
  inner_join(dat %>% filter(treatment != "unirradiated")) %>%
  mutate(ddCT = dCT - ctrl_dCT) %>%
  ungroup()

dat %>%
  anova_test(dCT~treatment*gene)

## plot data
svg("Results/Tissue_paper/cSCC_mRNA_irradiation.svg",  width=5, height=4)
ddCT %>% 
  ggbarplot(x = "gene", y = "ddCT", add = "mean_se", fill = "grey40") +
  geom_jitter(alpha = 0.5, position = position_dodge(0.6)) + 
  geom_hline(yintercept = 0, lty = 3) + 
  # theme_PhD(axis.text.size = 12) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  scale_y_continuous(guide = "axis_minor",limits = c(-3.59,5.56), breaks = seq(-4,6,2)) +
  ylab("log2 fold-change of gene expression")
dev.off()

## plot data
svg("Results/Tissue_paper/cSCC_CDH_ratio_irradiation.svg",  width=2, height=4)
EMT_markers %>% 
  ggbarplot(x = "treatment", y = "ratio", add = "mean_se", fill = "treatment") +
  geom_jitter(alpha = 0.5, position = position_dodge(0.6)) + 
  geom_hline(yintercept = 0, lty = 3) + 
  scale_fill_manual(values = c("white", "grey40")) +
  # theme_PhD(axis.text.size = 12) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  scale_y_continuous(guide = "axis_minor") +
  ylab("Ratio of relative CDH2/CDH1 expression")
dev.off()


