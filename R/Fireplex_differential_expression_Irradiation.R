## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(ggpubr)
library(rstatix)
library(agricolae)
library(EnvStats)
library(devtools)
library(data.table)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

#load data
url_file <- "https://raw.githubusercontent.com/MBender1992/PhD/Marc/Data/PhD_MB_FirePlex_chronic_irr_20190620.csv" 
dat <-  load_Fireplex_data_PhD(filename = url(url_file), threshold = 2.5)


###################
# Data inspection #
###################

# summary statistics 
dat_summary <- dat %>% 
  mutate(expression = expression + 1) %>%
  group_by(cell_line, Irradiation,miRNA) %>%
  summarize(geomean = geoMean(expression, na.rm =T), geosd = geoSD(expression))


# visualize data
dat %>% 
  ggboxplot(
    x = "cell_line",
    y = "expression",
    color = "Irradiation",
    facet.by = "miRNA",
    scales = "free", 
    palette = "jco"
  )



########################
# Checking assumptions #
########################

# check for outliers 
dat %>% 
  group_by(cell_line,Irradiation, miRNA) %>%
  identify_outliers(log_exp) 

# normality assumption
model <- lm(log_exp~cell_line*miRNA*Irradiation, data = dat)
ggqqplot(residuals(model))


# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

dat %>%
  group_by(miRNA, Irradiation,cell_line) %>%
  shapiro_test(log_exp)

ggqqplot(dat, "log_exp", ggtheme = theme_bw()) +
  facet_grid(cell_line + Irradiation ~ miRNA, labeller = "label_both")

# homogeneity of variance assumption
dat %>%
  levene_test(log_exp~cell_line*miRNA*Irradiation)

# homogeneity of variance within each miRNA
dat %>% 
  group_by(miRNA) %>%
  levene_test(log_exp~cell_line*Irradiation) %>% arrange(p) %>%print(n="all")


#################
#   Statistics  #
#################

# expression values are log normal distributed
# to account for heterogeneity of variances use Welch ANOVA

## =================================
## ANOVA 
# 3-way ANOVA
dat%>% 
  anova_test(log_exp ~ miRNA*cell_line*Irradiation, white.adjust = T)

# 2-way ANOVA
two_way_ANOVA <- dat%>% 
  group_by(miRNA) %>%
  anova_test(log_exp ~ cell_line*Irradiation, white.adjust = T) %>%
  adjust_pvalue(method = "fdr")

# significant interaction effect
two_way_interaction <- two_way_ANOVA %>%
  filter(Effect == "cell_line:Irradiation" & p.adj < 0.05)%>%
  pull(miRNA)
 
# significant Irradiation main effect
two_way_irradiation <- two_way_ANOVA %>%
  filter(Effect == "Irradiation" & p.adj < 0.05 & !miRNA %in% two_way_interaction) %>%
  pull(miRNA)



## =================================
## Post-hoc Tests 

## Interaction effect
# calculate fold changes
folds <-  dat_summary %>% 
  group_by(cell_line,  miRNA) %>% 
  summarize(FC = geomean[Irradiation == "KAUVIR"]/geomean[Irradiation == "control"]) %>% 
  mutate(log2FC = log2(FC)) %>%
  setNames(c("cell_line", "miRNA","FC", "log2FC"))

# t-test for the interaction effect
pwc_interaction <- dat %>%
  group_by(miRNA,cell_line) %>%
  pairwise_t_test(log_exp~Irradiation, pool.sd = F) %>% 
  filter(p < 0.05 & miRNA %in% two_way_interaction) 

#calculate differential expression for interaction
diff_expr_interaction <- pwc_interaction %>%
  inner_join(folds) %>%
  select(miRNA, cell_line, p.adj, p.adj, log2FC) %>% 
  filter(abs(log2FC) > log2(1.5))



## Irradiation effect

# post-hoc-test for the Irradiation effect
pwc_irradiation <- dat%>%
  group_by(miRNA,cell_line) %>%
  pairwise_t_test(log_exp~Irradiation, pool.sd = F) %>% 
  filter(p < 0.05 & miRNA %in% two_way_irradiation) 

#calculate differential expression for Irradiation
diff_expr_irradiation <- pwc_irradiation %>%
  inner_join(folds) %>%
  select(miRNA, cell_line, p.adj, p.adj, log2FC) %>% 
  filter(abs(log2FC) > log2(1.5))



##################################
#    Prepare data for plotting   #
##################################


# expand limits for each facet
dat_plot <- facet_limits(dat_summary, "geomean + geosd", "miRNA") %>% 
  mutate(cell_line = str_replace_all(cell_line, "_","-"))

# define limits of errorbar to only show upper errorbar
limits <- aes(ymax = ifelse(geomean>0,geomean + geosd,geomean/2),
              ymin = ifelse(geomean<0,geomean - geosd,geomean/2))




##################################
#    Plot interaction effect     #
##################################

# plot interaction effects of cell line and irradiation (irradiation effect differs depending on cell context)
png("Results/Fireplex_Irradiation_Interaction.png", units="in", width=7.5, height=6, res=600)


dat_plot %>% 
  filter(miRNA %in% diff_expr_interaction$miRNA) %>%
  ggplot(
    aes(
      x=cell_line,
      y=geomean,
      fill=Irradiation
      )
    ) + 
  geom_errorbar(
    limits,
    width=.2,
    position=position_dodge(0.5)
    ) +
  geom_bar(
    stat="identity",
    color="black",
    position=position_dodge(),
    width= 0.5
    ) +
  facet_wrap(~miRNA,scales="free")+ 
  theme_PhD(axis.text.size = 8) +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    legend.key.size = (unit(0.7,"line")),
    legend.position = "bottom"
  ) +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max)) + 
  scale_fill_manual(values = grey.colors(5,start=1,end=0.2))  +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("miRNA expression (a.u.)")

dev.off()



##################################
#    Plot Irradiation effect     #
##################################

# plot interaction effects of cell line and irradiation (irradiation effect differs depending on cell context)
png("Results/Fireplex_Irradiation_Main_Effect.png", units="in", width=10, height=6, res=600)


dat_plot %>% 
  filter(miRNA %in% diff_expr_irradiation$miRNA) %>%
  ggplot(
    aes(
      x=cell_line,
      y=geomean,
      fill=Irradiation
    )
  ) + 
  geom_errorbar(
    limits,
    width=.2,
    position=position_dodge(0.5)
  ) +
  geom_bar(
    stat="identity",
    color="black",
    position=position_dodge(),
    width= 0.5
  ) +
  facet_wrap(~miRNA,scales="free")+ 
  theme_PhD(axis.text.size = 8) +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    legend.key.size = (unit(0.7,"line")),
    legend.position = "bottom"
  ) +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max)) + 
  scale_fill_manual(values = grey.colors(5,start=1,end=0.2))  +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("miRNA expression (a.u.)")

dev.off()

