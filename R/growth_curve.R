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
dat <- readxl::read_excel("Data/20212905_growth_curve_transfection_antago_mir_205.xlsx")

df <- dat %>% 
  group_by(time, treatment) %>%
  summarize(mean = mean(cell_number), sd = sd(cell_number)) %>%
  mutate(lower = mean - sd, upper = mean + sd) %>%
  mutate(treatment = ifelse(treatment == "ctrl", "Kontrolle", "miR-205 Inhibitor"))



dat%>%
  anova_test(cell_number~treatment*time)



svg("Results/Growth_curve.svg", width=4, height=3)
df %>% 
  ggplot(aes(time, mean, color = treatment, shape = treatment)) + 
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.7, width = 3) + 
  geom_line() + 
  theme_PhD(axis.text.size = 12) +
  theme(
    strip.background = element_blank(),
    legend.key.size = (unit(0.7,"line")),
    legend.position = "bottom"
  ) +
  scale_color_jco() +
  scale_x_continuous(limits = c(0,172),breaks = seq(0, 168, 24),expand = c(0, 0), guide = "axis_minor") +
  scale_y_continuous(limits = c(0,650000), expand = c(0, 0), breaks = seq(0,600000,100000), guide = "axis_minor") +
  xlab("Zeit nach Aussaat (h)") + 
  ylab("Anzahl Zellen")
dev.off()


