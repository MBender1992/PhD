## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(ggpubr)
library(ggsci)
library(devtools)
library(ggh4x)
library(rstatix)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R") 

# load data 
dat <- read_csv("Data/EMT_marker_extracellular_v0.2_20220404.csv") %>%
  filter(!Messung %in% c("20220120")) %>%
  mutate(Name = factor(Name, levels = c("negative_control", "transfection_control", "anti_miR_205", "TGF_beta_10ng", "TGF_beta_10ng_48h", "TGF_beta_130ng", "anti_miR_205+TGF_beta", "fibroblasts"))) %>%
  select(-Messung)
colnames(dat) <- c("ID", "Name", "N-Cadherin", "E-Cadherin")

svg("Results/EMT_markers.svg",  width=3, height=5)
dat %>% gather(-c("ID", "Name"), key = "marker", value ="expression") %>%
ggplot(aes(marker, expression, fill = marker)) +
  geom_bar(stat = "summary", fun = "mean",
           color = "black", position = position_dodge(0.6),
           mapping = aes(x = marker,y = expression),
           width=0.6
  )  +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
                position = position_dodge(0.6), width = 0.25, size = 0.6) +
  facet_wrap(~Name, scales = "free", nrow = 4) +
  # geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(0.5)) +
  # geom_jitter(alpha = 0.7, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)) +
  theme_PhD(axis.text.size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) + 
  scale_y_continuous(expand = c(0,0), guide = "axis_minor") +
  scale_fill_jco(alpha = 0.7) +
  #scale_fill_brewer(palette = "Greys")+
  ylab("Proteinexpression (a.u.)")
dev.off()


# display ratio of E to N-Cadherin
dat_ratio <- dat %>% mutate(ratio = log2(`N-Cadherin`/`E-Cadherin`)) 

# Dunnett Test to compare means vs control
anova_test(dat_ratio, ratio~Name, white.adjust = TRUE)
DescTools::DunnettTest(dat_ratio$ratio, dat_ratio$Name, control = "negative_control", conf.level = 0.95)

svg("Results/EMT_markers_ratio.svg",  width=3, height=3)
dat_ratio %>%
  # mutate(Name = fct_reorder(Name, ratio)) %>%
  ggplot(aes(Name, ratio, fill =Name)) + 
  geom_boxplot() +
  # geom_bar(stat = "summary", fun = "mean",
  #          color = "black", position = position_dodge(0.6),
  #          mapping = aes(x = Name,y = ratio),
  #          width=0.6)   +
  # geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
  #               position = position_dodge(0.6), width = 0.25, size = 0.6) +
  theme_PhD(axis.text.size = 12) +
  theme(
    axis.title.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    strip.background = element_blank(), 
    legend.position = "none",
    legend.title = element_blank()
  )  + 
  ylab("N-Cadherin/E-Cadherin (log2)")+
  # geom_hline(yintercept = 0, lty = 1) +
  scale_y_continuous(guide = "axis_minor", limits = c(-1.7, 2.7), breaks = c(seq(-1.5,2.5, 0.5))) +
  scale_fill_brewer(palette = "Blues")
dev.off()
 


