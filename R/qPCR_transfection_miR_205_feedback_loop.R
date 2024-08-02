## <<<<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(devtools)
library(EnvStats)
library(readxl)
library(ggh4x)
library(data.table)
library(rstatix)
library(ggsci)


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
  ungroup() %>%
  rename(Zeitpunkt = time)

# statistics
ddCT %>% group_by(gene_name, Zeitpunkt) %>% 
  t_test(ddCT~ 1, mu = 0) %>% 
  adjust_pvalue(method = "fdr") %>%
  # filter(p <= 0.05) %>% 
  arrange(gene_name)




# prepare data for plotting
dat_summary <- ddCT %>% 
  group_by(Zeitpunkt, gene_name) %>%
  summarize(mean = mean(ddCT, na.rm =T), sd = sd(ddCT)) %>%
  mutate(Zeitpunkt = factor(Zeitpunkt, levels = c("24h", "72h"))) %>%
  mutate(ymin = mean - sd, ymax = mean + sd)

# dat_plot <- facet_limits(dat_summary, "mean + sd", "gene_name") 
# dat_points <- ddCT %>% rename(mean = ddCT)

# define limits of errorbar to only show upper errorbar
# limits <- aes(ymax = ifelse(mean>0,mean + sd,mean),
#               ymin = ifelse(mean<0,mean - sd,mean))  

svg("Results/Antago_miR-205_feedback.svg", width=2.5, height=4)
dat_summary %>% 
  ggplot(aes(gene_name, mean, fill = Zeitpunkt)) +
  geom_bar(stat="identity",color="black",  position=position_dodge(0.5),  width= 0.5)+
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width=.2,  position=position_dodge(0.5)) +
  # geom_point(data = dat_points, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5), alpha = 0.5) +
  geom_hline(yintercept= 0, lty = 2) +
  theme_PhD(axis.text.size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.background = element_blank(),  legend.position = "bottom") + 
  scale_fill_jco(alpha = 0.7) +
  scale_y_continuous(limits = c(-4.5, 2.5), breaks = seq(-4,2,1), guide = "axis_minor") +
  ylab("ddCT")
dev.off()


