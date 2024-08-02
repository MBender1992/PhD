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
dat <- read_csv("Data/Results_scratch_assay_transfection.csv") %>%
  mutate(wound_closure = wound_closure*100)


pixel_to_µm <- 1772/421.94
pixel2_to_µm2 <- 3605314/204391
dat <- dat %>% mutate(Area_µm = Area_pixels/pixel2_to_µm2)

# dat <- dat %>% group_by(Messung, time, treatment) %>% summarize(wound_closure = mean(wound_closure))
dat <- dat %>% filter(!is.na(Area_µm))
dat_summary <- Rmisc::summarySE(dat, measurevar="Area_µm", groupvars=c("time","treatment"))

svg("Results/scratch_assay_miR205_Woundclosure.svg", width=4, height=3)
dat_summary %>% 
  ggplot(aes(time, wound_closure, color = treatment)) + 
  geom_smooth(method = "lm", se = F, lty = 1) +
  # geom_line() +
  geom_errorbar(aes(ymin = wound_closure - sd, ymax = wound_closure + sd), width = 0.4) +
  geom_point(size = 2) +
  theme_PhD(axis.text.size = 12) +
  theme(
    axis.title.x = element_text(face = "bold", size = 8),
    strip.background = element_blank(),
    legend.key.size = (unit(0.7,"line")),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_jco() +
  xlab("Time (h)") +
  ylab("Wound closure (%)") + 
  scale_x_continuous(breaks = seq(0,24,4), guide = "axis_minor") +
  scale_y_continuous(expand = c(0, 0), guide = "axis_minor", breaks = seq(0, 100, 20))
dev.off()  


svg("Results/scratch_assay_miR205_Area.svg", width=4, height=3)
dat_summary %>% 
  ggplot(aes(time, Area_µm, color = treatment)) + 
  geom_smooth(method = "lm", se = F, lty = 1) +
  # geom_line() +
  geom_errorbar(aes(ymin = Area_µm - sd, ymax = Area_µm + sd), width = 0.4) +
  geom_point(size = 2) +
  theme_PhD(axis.text.size = 12) +
  theme(
    axis.title.x = element_text(face = "bold", size = 8),
    strip.background = element_blank(),
    legend.key.size = (unit(0.7,"line")),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_jco() +
  xlab("Time (h)") +
  ylab(expression(bold(Fläche~("\u03bc"~m^2))))  +
  scale_x_continuous(breaks = seq(0,24,4), guide = "axis_minor") +
  scale_y_continuous(expand = c(0, 0), guide = "axis_minor", limits = c(-40000, 280000), breaks = seq(0, 250000, 50000))
dev.off()  

# RMD Report?
# Primer