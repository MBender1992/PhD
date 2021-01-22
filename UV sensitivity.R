library(tidyverse)
library(ggpubr)
library(rstatix)
library(export)
library(ggrepel) 
library(car)
library(devtools)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  


#############################
#                           #
# 1. Load and process data  #
#                           #
#############################

#load data
url_file <- "https://raw.githubusercontent.com/MBender1992/PhD/Marc/Data/200619_chronic_irr_normalized.csv" 
dat <-  load_Fireplex_data_PhD(filename = url(url_file), threshold = 2.5)

# summarize data
dat_summarized <- dat %>% 
  group_by(cell_line, Irradiation, miRNA) %>% 
  summarize(mean = mean(expression)) %>%
  ungroup() %>% 
  group_by(cell_line, miRNA) %>%
  spread(Irradiation, mean) %>%
  mutate(FC = ifelse(is.infinite(KAUVIR/control), NA, KAUVIR/control)) %>%
  mutate(logFC = ifelse(is.infinite(log2(FC)), NA, log2(FC))) %>%
  mutate(abslogFC = abs(logFC)) %>% 
  ungroup() %>% 
  mutate(miRNA = str_replace_all(.$miRNA, "hsa_","")) %>%
  mutate(miRNA = str_replace_all(.$miRNA, "_","-")) 
  


#################################
#                               #
# 2. Explorative data analysis  #
#                               #
#################################


########################
# Checking assumptions #
########################

# check for outliers 
dat_summarized %>% 
  group_by(cell_line) %>%
  identify_outliers(abslogFC)

# normality assumption
model <- lm(abslogFC~cell_line, data = dat_summarized)
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# shapiro test on each group
dat_summarized %>%
  group_by(cell_line) %>%
  shapiro_test(abslogFC)

# qqplot for each group
ggqqplot(dat_summarized, "abslogFC", ggtheme = theme_bw()) +
  facet_wrap(~ cell_line)

# homogeneity of variance assumption
dat_summarized %>% 
  levene_test(abslogFC~cell_line)



#### assumption of normality and homogeneity of variance violated

# ..................................................................................................................................
# dat_summarized %>% welch_anova_test(abslogFC~cell_line)
# pvals <- dat_summarized %>% 
#   games_howell_test(abslogFC~cell_line)

dat_summarized %>% 
  wilcox_test(abslogFC~cell_line, p.adjust.method = "bonferroni")


################
#              #
# 3. Plotting  #
#              #
################


# calculate and plot UV sensitivity as mean of absolute log Fold-changes 
# define position for jitter and text labels
pos <- position_jitter(width = 0.05, seed = 2)

# define labels
labels <- dat_summarized %>% 
  group_by(cell_line) %>%
  top_n(2, abslogFC) %>% 
  mutate(label = miRNA)

# merge data with labels
dat_labelled <- left_join(dat_summarized,labels) %>%
  mutate(label = ifelse(is.na(label), "", label))

# open png device
png("test.png", units="in", width=5.5, height=6, res=600)

# plot
dat_labelled %>% 
  ggplot(aes(cell_line,abslogFC)) +
  stat_boxplot(
    geom ='errorbar',
    size=0.5,
    width = 0.2
  ) +
  geom_boxplot(
    color = "black",
    mapping = aes(x = cell_line,y = abslogFC),
    outlier.shape =NA,
    width=0.4
  ) +
  geom_jitter(
    data =  filter(dat_labelled, label == ""),
    mapping = aes(x = cell_line,y = abslogFC),
    width=0.1,
    alpha=0.3
  ) +
  geom_jitter(
    data =  filter(dat_labelled, label != ""),
    mapping = aes(x = cell_line,y = abslogFC, fill = "red"),
    shape = 21, 
    colour = "black",
    position = pos
  ) + 
  geom_text_repel(
    data =  filter(dat_labelled, label != ""),
    mapping = aes(x = cell_line, y = abslogFC, label = label),
    size = 3,
    position = pos
  ) + 
  # stat_summary(fun = mean, geom = "crossbar", width = 0.4, lwd = 0.25) +
  # stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) +
  theme_PhD(Legend=F) +
  theme(axis.text.x = element_text(angle=90, vjust= 0.4, hjust= 1)) +
  ylab("absolute log-FC") +
  scale_y_continuous(limits = c(0,3.3), breaks = seq(0,3,0.5))

dev.off()

# save file with transparency as power point
# graph2ppt(file="UV_sensitivity.pptx", width=3.5, height=5)
  






