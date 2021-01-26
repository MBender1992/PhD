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
dat_UV_sens <- folds %>% 
  ungroup()%>% 
  mutate(abslog2FC= abs(log2FC))
  

#################################
#                               #
# 2. Explorative data analysis  #
#                               #
#################################


########################
# Checking assumptions #
########################

# check data
dat_UV_sens %>% 
  ggplot(aes(cell_line, abslog2FC)) +
  geom_boxplot()

# check for outliers 
dat_UV_sens %>% 
  group_by(cell_line) %>%
  identify_outliers(abslog2FC)

# normality assumption
model <- lm(abslog2FC~cell_line, data = dat_UV_sens)
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# shapiro test on each group
dat_UV_sens %>%
  group_by(cell_line) %>%
  shapiro_test(abslog2FC)

# qqplot for each group
ggqqplot(dat_UV_sens, "abslog2FC", ggtheme = theme_bw()) +
  facet_wrap(~ cell_line)

# homogeneity of variance assumption
dat_UV_sens %>% 
  levene_test(abslog2FC~cell_line)



#### assumption of normality and homogeneity of variance violated

# ..................................................................................................................................
# dat_UV_sens %>% welch_anova_test(abslog2FC~cell_line)
# pvals <- dat_UV_sens %>% 
#   games_howell_test(abslog2FC~cell_line)

# Kruskal Wallis test to test if there are any differences beteen cell lines
dat_UV_sens %>% 
  kruskal_test(abslog2FC~cell_line)

# Wilcoxon-Mann-Whitney U test to compare groups pairwise
dat_UV_sens %>% 
  wilcox_test(abslog2FC~cell_line, p.adjust.method = "fdr")


################
#              #
# 3. Plotting  #
#              #
################


# calculate and plot UV sensitivity as mean of absolute log Fold-changes 
# define position for jitter and text labels
pos <- position_jitter(width = 0.05, seed = 2)

# define labels
labels <- dat_UV_sens %>% 
  group_by(cell_line) %>%
  top_n(3, abslog2FC) %>% 
  mutate(label = miRNA)

# merge data with labels
dat_labelled <- left_join(dat_UV_sens,labels) %>%
  mutate(label = ifelse(is.na(label), "", label))

# open png device
png("test.png", units="in", width=5.5, height=6, res=600)

# plot
dat_labelled %>% 
  ggplot(aes(cell_line,abslog2FC)) +
  stat_boxplot(
    geom ='errorbar',
    size=0.5,
    width = 0.2
  ) +
  geom_boxplot(
    color = "black",
    mapping = aes(x = cell_line,y = abslog2FC),
    outlier.shape =NA,
    width=0.4
  ) +
  geom_jitter(
    data =  filter(dat_labelled, label == ""),
    mapping = aes(x = cell_line,y = abslog2FC),
    width=0.1,
    alpha=0.3
  ) +
  geom_jitter(
    data =  filter(dat_labelled, label != ""),
    mapping = aes(x = cell_line,y = abslog2FC, fill = "red"),
    shape = 21, 
    colour = "black",
    position = pos
  ) + 
  geom_text_repel(
    data =  filter(dat_labelled, label != ""),
    mapping = aes(x = cell_line, y = abslog2FC, label = label),
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
  






