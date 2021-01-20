##Pr#ambel
setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/base_scripts")

#loading custom functions
source("R_functions.R")
# source("R_functions_PhD.R")

library("tidyverse")
library("ggpubr")
library("rstatix")
library("EnvStats")
library("data.table")
library("agricolae")

setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/PhD/Daten")

#load data
filename <- "200619_chronic_irr_normalized.csv"
dat <-  load_Fireplex_data_PhD(threshold = 2.5)


###################
# Data inspection #
###################

# summary statistics 
dat_summary <- dat %>% 
  mutate(expression = expression +0.1) %>%
  group_by(cell_line, Irradiation,miRNA) %>% 
  summarize(geomean = geoMean(expression,na.rm=T),
            geosd = geoSD(expression,na.rm=T))

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

dat%>%
  group_by(miRNA, Irradiation,cell_line) %>%
  shapiro_test(log_exp)

ggqqplot(dat, "log_exp", ggtheme = theme_bw()) +
  facet_grid(cell_line + Irradiation ~ miRNA, labeller = "label_both")

# homogeneity of variance assumption
dat%>%
  levene_test(log_exp~cell_line*miRNA*Irradiation)

# homogeneity of variance within each miRNA
dat%>% 
  group_by(miRNA) %>%
  levene_test(log_exp~cell_line*Irradiation) %>% arrange(p) %>%print(n="all")


#################
#   Statistics  #
#################

# expression values are log normal distributed
# to account for heterogeneity of variances use Welch ANOVA


## ANOVA .........................................................................................................
# 3-way ANOVA
dat%>% 
  anova_test(log_exp ~ miRNA*cell_line*Irradiation, white.adjust = T)

# 2-way ANOVA
two_way_ANOVA <- dat%>% 
  group_by(miRNA) %>%
  anova_test(log_exp ~ cell_line*Irradiation, white.adjust = T) %>%
  adjust_pvalue(method = "fdr")


# significant interaction effect
two_way_ANOVA_interaction <- two_way_ANOVA %>%
  filter(Effect == "cell_line:Irradiation" & p.adj < 0.05) 
 
# significant Irradiation main effect
two_way_ANOVA_Irradiation <- two_way_ANOVA %>%
  filter(Effect == "Irradiation" & p.adj < 0.05 & !miRNA %in% two_way_ANOVA_interaction$miRNA) 


## Post-hoc Tests ................................................................................................
# post-hoc-test for the interaction effect
pwc_interaction <- dat%>%
  group_by(miRNA,cell_line) %>%
  pairwise_t_test(log_exp~Irradiation, pool.sd = F) %>% 
  filter(p.adj < 0.05 & miRNA %in% two_way_ANOVA_interaction$miRNA) 


# post-hoc-test for the Irradiation effect
pwc_Irradiation <- dat%>%
  group_by(miRNA,cell_line) %>%
  pairwise_t_test(log_exp~Irradiation, pool.sd = F) %>% 
  filter(p.adj < 0.05 & miRNA %in% two_way_ANOVA_Irradiation$miRNA) 




## Differential expression..........................................................................................
# calculate fold changes
folds <- dat_summary %>% 
  group_by(cell_line,  miRNA) %>% 
  summarize(FC = geomean[Irradiation == "KAUVIR"]/geomean[Irradiation == "control"]) %>%
  mutate(logFC = log2(FC))

#calculate differential expression for interaction
diff_expr_interaction <- pwc_interaction %>%
  inner_join(folds) %>%
  select(miRNA, cell_line, p.adj, p.adj, FC,logFC) %>% 
  filter(abs(logFC) > log2(1.5))

#calculate differential expression for Irradiation
diff_expr_Irradiation <- pwc_Irradiation %>%
  inner_join(folds) %>%
  select(miRNA, cell_line, p.adj, p.adj, FC,logFC) %>% 
  filter(abs(logFC) > log2(1.5))




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
png("test.png", units="in", width=7.5, height=6, res=600)


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
png("test.png", units="in", width=10, height=6, res=600)


dat_plot %>% 
  filter(miRNA %in% diff_expr_Irradiation$miRNA) %>%
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

