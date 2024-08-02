## <<<<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(devtools)
library(readxl)
library(rstatix)
library(ggh4x)
library(ggsci)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

dat <- read.csv("Data/Invasion_Assay_Boyden_Chamber_20220823.csv")

dat$Name <- with(dat, paste0(upper_chamber, lower_chamber,treatment))
dat$Name <- factor(dat$Name, levels = c("++ctrl", "-+ctrl", "-+TGFb"))

dat_ratio <- dat %>% 
  group_by(Name, replicate) %>%
  summarize(ratio = mean[membrane == "bottom"]/mean[membrane == "all"]) %>%
  ungroup()

dat$membrane <- factor(dat$membrane, levels = c("all", "bottom"), labels = c("Gesamtanzahl Zellen", "Eingedrungene Zellen"))

svg("Results/Invasion_Assay_raw.svg", width=6, height=2.5)
dat %>% 
  ggplot(aes(membrane, mean, fill = membrane)) + 
  facet_wrap(~Name, nrow = 1, scales = "free") +
  geom_bar(stat = "summary", fun = "mean",
           color = "black", position = position_dodge(0.7), 
           width = 0.6) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
                position = position_dodge(0.6), width = 0.25, size = 0.6) +
  theme_PhD(axis.text.size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )+
  scale_y_continuous(expand = c(0,0), limits = c(0,1200), breaks = seq(0,1200,200), guide = "axis_minor") +
  scale_fill_jco(alpha = 0.7) +
  ylab("Anzahl Zellen pro Sichtfeld")
dev.off()

svg("Results/Invasion_Assay_ratio.svg",  width=4, height=3)
dat_ratio %>% 
  ggplot(aes(Name, ratio, fill = Name))+
  geom_bar(stat = "summary", fun = "mean",
           color = "black", position = position_dodge(0.7), 
           width = 0.6) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
                position = position_dodge(0.6), width = 0.25, size = 0.6)+
  theme_PhD(axis.text.size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1,0.2), guide = "axis_minor", labels = scales::percent) +
  scale_fill_brewer(palette = "Greys") +
  ylab("Eingedrungene  Zellen")
dev.off()


test <- dat %>% group_by(Name, membrane) %>% summarize(mean = mean(mean)) %>%
  spread(membrane, mean) %>%
  setNames(c("Name", "n", "invaded")) %>%
  mutate(non_migrated = n-invaded) 

# global chi-square test
chisq_test(test[,3:4])

# groupwise comparison
# ++ctrl vs -+ctrl
chisq_test(test[-3,3:4])
p1<-prop.test(test$invaded[-3], test$n[-3])

# ++ctrl vs TGF
chisq_test(test[-2,3:4])
p2<-prop.test(test$invaded[-2], test$n[-2])

# -+ctrl vs TGF
chisq_test(test[-1,3:4])
p3<-prop.test(test$invaded[-1], test$n[-1])

p<-c(p1$p.value,p2$p.value,p3$p.value) 

p.adjust(p, method = "fdr", n = length(p))



