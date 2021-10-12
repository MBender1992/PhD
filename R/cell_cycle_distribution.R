## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(ggpubr)
library(ggsci)
library(devtools)
library(ggh4x)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R") 

# load data 
dat <- readxl::read_excel("Data/20210601_cell_cycle_transfection_antago_mir_205.xlsx") %>%
  gather(key = "phase", value = "percent", -c(ID, treatment, time)) %>%
  mutate(time = factor(time, levels = c("24 h", "48 h", "72 h", "96 h", "168 h"))) %>%
  mutate(phase = factor(phase, levels = c("G2", "S", "G1"))) %>% 
  mutate(treatment = factor(treatment, levels = c("control", "antagomir-205")))

png("Results/cell_cycle_distribution.png", units="in", width=6, height=3, res=300)
dat %>%
  ggbarplot(x = "time", y = "percent", fill = "phase", alpha = 0.7,
            add = "mean_se", facet.by = "treatment", scales = "free") + 
  theme_PhD(axis.text.size = 8) +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    legend.key.size = (unit(0.7,"line")),
    legend.position = "bottom"
  ) +
  scale_fill_jco() +
  scale_y_continuous(expand = c(0, 0), guide = "axis_minor", breaks = seq(0,100,20)) +
  ylab("% of gated cells")
dev.off()
