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

dat_ELISA <- read.csv("Data/ELISA_20230307.csv") 

fit1 <- lm(Konzentration ~ OD, data = dat_ELISA[dat_ELISA$Messung == 1,])
fit2 <- lm(Konzentration ~ OD, data = dat_ELISA[dat_ELISA$Messung == 2,])

dat_ELISA <- dat_ELISA %>%
  filter(Type != "Standard") %>%
  mutate(Konzentration = ifelse(Messung == 1, (fit1$coefficients[1]+fit1$coefficients[2]*OD)*30/1000, (fit2$coefficients[1]+fit2$coefficients[2]*OD)*30/1000)) %>%
  mutate(Konzentration = ifelse(Konzentration < 0, 0, Konzentration)) %>%
  mutate(Bestrahlung = factor(Bestrahlung, level = c("Kontrolle", "Irr"))) 
# multipliziert mit 30 für den Verdünnungsfaktor
# geteilt durch 1000 um die Konzentration in ng/mL zu erhalten


dat_summary <- dat_ELISA %>% filter(Messung == 2) %>% 
  group_by(Zeitpunkt, Zelltyp, Bestrahlung) %>%
  summarize(mean = mean(Konzentration, na.rm =T), sd = sd(Konzentration)) %>%
  mutate(Zeitpunkt = factor(Zeitpunkt, levels = c("4h", "72h"))) %>%
  mutate(ymin = mean - sd, ymax = mean + sd)

svg("Results/Tgfb_ELISA.svg")
dat_summary %>%
  ggplot(aes(Zeitpunkt, mean, fill = Bestrahlung)) + facet_wrap(~Zelltyp, scales = "free") + 
  geom_bar(stat="identity",color="black",  position=position_dodge(0.5),  width= 0.5)+
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width=.2,  position=position_dodge(0.5)) +
  theme_PhD(axis.text.size = 10) +
  theme(
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),
    legend.position = "none"
  ) + scale_y_continuous(expand = c(0,0), limits = c(0,2), breaks = seq(0,2,0.5)) +scale_fill_jco(alpha = 0.7)
dev.off()

dat_ELISA %>% filter(Messung == 2) %>% group_by(Zelltyp, Zeitpunkt) %>%
  tukey_hsd(Konzentration ~ Bestrahlung)


