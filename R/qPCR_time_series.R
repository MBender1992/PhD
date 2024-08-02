##Pr#ambel
setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/base_scripts")

#loading custom functions
source("R_functions.R")
# source("R_functions_PhD.R")

library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)
library(EnvStats)

setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/PhD/Daten")


# load data and calculate mean of technical replicates
dat_qPCR <- read_xlsx("Data/qPCR_time_series_miR_205.xlsx", sheet = "Results")  %>%
  select(-Messung_1) %>%
  mutate(
    treatment = str_replace_all(.$Name,"^\\d+_","") %>%
      str_replace_all(".?\\d{1}x","") %>%
      str_replace_all("[1-9]$", "") %>%
      str_replace_all("K0", "con")
  ) %>% 
  mutate(
    dose = ifelse(str_detect(.$Name, "1x"), "1x250 J/m²", ""),
    dose = ifelse(str_detect(.$Name, "5x"), "5x250 J/m²", dose),
    dose = ifelse(str_detect(.$Name, "8x"), "8x250 J/m²", dose)
   ) %>%
  mutate(dose = factor(dose, levels = c("1x250 J/m²","5x250 J/m²","8x250 J/m²"))) 



  

# calculate geoMean of HK genes for each separate petri dish
dat_HK <- dat_qPCR %>% 
  group_by(Name, cell_line) %>%
  filter(gene_type == "HK") %>%
  mutate(geomean_HK = geoMean(mean_Ct, na.rm=T)) %>%
  distinct(geomean_HK)

# calculate dCT values and merge dCT values of technical replicates 
dat_dCT <- inner_join(dat_qPCR, dat_HK) %>%
  ungroup() %>%
  filter(gene_type != "HK") %>%
  select(-gene_type) %>%
  mutate(dCT = geomean_HK - mean_Ct)


## calculate ddCT
ddCT <- dat_dCT %>% 
  group_by(gene_name,cell_line, dose) %>% 
  filter(treatment == "con") %>%
  summarize(ctrl_dCT = mean(dCT, na.rm = T)) %>%
  inner_join(dat_dCT %>% filter(treatment == "irr")) %>%
  mutate(ddCT = dCT - ctrl_dCT) %>% 
  ungroup()



######################
# Statistics for ddCT #
######################

# assumptions ddCT
model <- lm(ddCT~dose, data = ddCT) # add cell_line later
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# one-sided t-test vs null
ddCT %>% group_by(dose) %>%
  t_test(ddCT~ 1, mu = 0)

ddCT %>% 
  anova_test(ddCT~ dose)

nested_aov <- function(dat, x, y){
  m.full <- aov(dat[[y]]~dat[[x]])
  m.red  <- aov(dat[[y]]~0)
  anova(m.red, m.full)
}

# nested ANOVA 
nested_aov(ddCT, x = "dose", y = "ddCT")$`Pr(>F)`[2]

# Plotting
dat_dCT  %>%
  mutate(rel.expression = 2^dCT) %>%
  ggbarplot(
    x="dose",
    y = "rel.expression",
    fill = "treatment",
    position = position_dodge(0.6),
    add = "mean_sd",
    error.plot = "upper_errorbar",
    width = 0.6,
    scales = "free"
  ) +
  theme_PhD(axis.text.size = 8) +
  scale_fill_manual(values = grey.colors(5,start=1,end=0.2))  +
  ylab("relative expression") +
  scale_y_continuous(expand = c(0, 0)) 


# Plot log2 fold-changes
svg("Results/miR205_time_series.svg",  width=4, height=5)
ddCT  %>%
  filter(!is.na(ddCT)) %>%
  ggbarplot(
    x="dose",
    y = "ddCT",
    fill = "dose",
    add = "mean_sd",
    width = 0.5,
    scales = "free",
    error.plot = "errorbar",
  ) +
  theme_PhD(axis.text.size = 12, 
            Legend = F) +
  scale_fill_brewer("Blues")  +
  ylab("Log2-FC") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(limits = c(-1,0.2), breaks = seq(-1,0.2, 0.2))
  
dev.off()



# statistics on ddCT, and show results in log space


# dCT[trt] = ref[trt] - goi[trt]
# dCT[ctrl] = ref[ctrl] - goi[trl]
# ddCT = dCT[ctrl] - dCT[trt]

# anders als 2008 beschrieben (ist intuitiver da negative dCT Werte mit weniger Expression und positive dCT mit einer 
# h?heren Expression korrelieren)