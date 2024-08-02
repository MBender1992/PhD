## <<<<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(devtools)
library(EnvStats)
library(readxl)
library(rstatix)
library(ggh4x)


# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# load and transform data
dat_qPCR <- read_excel("Data/miR-205_Targets_240821.xlsx") %>% 
  rename(Name = Probenname) %>%
  mutate(time = str_extract(.$Name, "\\d+h")) %>%
  mutate(treatment = str_replace_all(.$Name, "_.+$","")) 


dat_HK <- dat_qPCR %>% 
  group_by(Name, Messung) %>%
  filter(gene_type == "HK") %>%
  mutate(geomean_HK = geoMean(mean2_ct, na.rm=T))  %>%
  distinct(geomean_HK) 

dat_dCT <- inner_join(dat_qPCR, dat_HK, by = c("Name", "Messung")) %>%
  ungroup() %>%
  filter(gene_type != "HK") %>%
  select(-gene_type) %>%
  mutate(dCT = geomean_HK - mean2_ct)

ddCT <- dat_dCT %>% 
  group_by(gene_name, time, Messung) %>%
  filter(treatment == "ctrl") %>%
  summarize(ctrl_dCT = mean(dCT, na.rm = T)) %>%
  inner_join(dat_dCT %>% filter(treatment != "ctrl")) %>%
  mutate(ddCT = dCT - ctrl_dCT) %>%
  ungroup()

# plot
dat_plot <- ddCT %>%
  mutate(time = factor(time, levels = c("24h", "72h", "144h")))


nested_aov <- function(dat, x, y){
  m.full <- aov(dat[[y]]~dat[[x]])
  m.red  <- aov(dat[[y]]~0)
  anova(m.red, m.full)
}

# nested ANOVA 
dat_24 <- dat_plot %>% filter(time == "24h")
dat_72 <- dat_plot %>% filter(time == "72h") 
dat_144 <- dat_plot %>% filter(time == "144h") 

p24 <- nested_aov(dat_24, x = "gene_name", y = "ddCT")$`Pr(>F)`[2]
p72 <- nested_aov(dat_72, x = "gene_name", y = "ddCT")$`Pr(>F)`[2]
p144 <- nested_aov(dat_144, x = "gene_name", y = "ddCT")$`Pr(>F)`[2]

p <- c(p24, p72, p144)

p.adjust(p, method = "fdr", n = length(p))


svg("Results/Antago_miR-205_targets.svg", width=7, height=5)
dat_plot %>% 
  filter(!is.na(ddCT)) %>% 
  ggplot(aes(gene_name, ddCT, fill = gene_name)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(0.1), alpha = 0.7) +
  facet_wrap(~time) +
  geom_hline(yintercept= 0, lty = 2) +
  theme_PhD(axis.text.size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),
    legend.position = "none"
  ) + 
  scale_y_continuous(limits = c(-2.2, 4), breaks = seq(-2,4, 1),  guide = "axis_minor") +
  scale_fill_brewer(palette = "Blues")
 dev.off()

 dat_plot %>% filter(gene_name == "MMP13")

 

 