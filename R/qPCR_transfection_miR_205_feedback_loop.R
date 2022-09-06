## <<<<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(devtools)
library(EnvStats)
library(readxl)
library(ggh4x)


# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# load and transform data
dat_qPCR <- read.csv("Data/PhD_MB_qPCR_transfection_mir_205_feedback_loop_20220602.csv") %>% 
  rename(Name = Probenname) %>%
  mutate(time = str_extract(.$Name, "\\d+h")) %>%
  mutate(treatment = str_replace_all(.$Name, "_.+$","")) 


dat_HK <- dat_qPCR %>% 
  group_by(Name) %>%
  filter(gene_type == "HK") %>%
  mutate(geomean_HK = geoMean(Ct, na.rm=T))  %>%
  distinct(geomean_HK) 

dat_dCT <- inner_join(dat_qPCR, dat_HK, by = c("Name")) %>%
  ungroup() %>%
  filter(gene_type != "HK") %>%
  select(-gene_type) %>%
  mutate(dCT = geomean_HK - Ct)

ddCT <- dat_dCT %>% 
  group_by(gene_name, time) %>%
  filter(treatment == "ctrl") %>%
  summarize(ctrl_dCT = mean(dCT, na.rm = T)) %>%
  inner_join(dat_dCT %>% filter(treatment != "ctrl")) %>%
  mutate(ddCT = dCT - ctrl_dCT) %>%
  ungroup()



# plot
dat_plot <- ddCT %>%
  mutate(time = factor(time, levels = c("24h", "72h")))

png("Results/Antago_miR-205_targets.png", units="in", width=7, height=5, res=600)
dat_plot %>% 
  filter(!is.na(ddCT)) %>% 
  ggplot(aes(gene_name, ddCT, fill = gene_name)) +
  # geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
  #               position = position_dodge(0.48), width = 0.25, size = 0.6) + # add only errorbars to one side of the bars
  # geom_bar(stat = "summary", fun = "mean",
  #          color = "black", position = position_dodge(0.48),
  #          mapping = aes(x = gene_name,y = ddCT),
  #          width=0.4
  # )  +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(0.1), alpha = 0.7) +
  facet_wrap(~time) +
  geom_hline(yintercept= 0, lty = 2) +
  theme_PhD(axis.text.size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),
    legend.position = "none"
  ) + 
  # scale_y_continuous(limits = c(-2.2, 4), breaks = seq(-2,4, 0.5)) +
  scale_fill_brewer(palette = "Blues")
dev.off()
