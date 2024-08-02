## <<<<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(devtools)
library(ggsci)
library(ggh4x)
library(rstatix)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

dat_apo <- read_csv("Data/141021_apoptosis_SCC12_transfection.csv") %>%
  mutate(transfection = factor(transfection, levels = c("ctrl", "anti-miR-205")),
         staurosporin = factor(staurosporin, levels = c("ctrl", "0.1 microM", "1 microM"), labels = c("0 \U00B5M", "0.1 \U00B5M", "1 \U00B5M"))) %>%
  select(-living_cells) %>%
  gather("apoptosis", "value", -c(ID, transfection, staurosporin)) %>%
  mutate(apoptosis = factor(apoptosis, levels = c("early_apoptosis", "late_apoptosis_necrosis"), labels = c("Frühapoptose", "Spätapoptose/Nekrose"))) %>%
  mutate(transfection = ifelse(transfection == "ctrl", "Kontrolle", "miR-205-5p Inhibitor")) 


# dat_apo %>%
#   group_by(staurosporin, apoptosis) %>%
#   t_test(value ~ transfection)
# 
# dat_apo %>% filter(apoptosis == "Early apoptosis") %>% 
#   t_test(value ~ staurosporin)
# 
# dat_apo %>%
#   anova_test(value~transfection*staurosporin*apoptosis) %>%
#   adjust_pvalue(method = "fdr")

svg("Results/apoptosis.svg",  width=8, height=3)
dat_apo %>%
  ggplot(aes(staurosporin, value , fill = transfection)) +
  geom_bar(stat = "summary", fun = "mean",
           color = "black", position = position_dodge(0.7),
           mapping = aes(x = staurosporin, y = value),
           width=0.6) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
                position = position_dodge(0.7), width = 0.25, size = 0.6) + # add only errorbars to one side of the bars
  facet_wrap(~ apoptosis, scales = "free")  +
  # geom_point(position = position_dodge(0.7)) +
  theme_PhD(axis.text.size = 12) +
  theme(strip.background = element_blank(), legend.position = "right") + 
  xlab("Staurosporin") +
  ylab("Anteil apoptotischer Zellen (%)") +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0,16, 2), expand = c(0, 0),guide = "axis_minor")+
  scale_fill_jco(alpha = 0.7)
dev.off()