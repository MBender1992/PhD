## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(ggpubr)
library(devtools)
library(rstatix)
library(EnvStats)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

#load data
url_file <- "https://raw.githubusercontent.com/MBender1992/PhD/Marc/Data/PhD_MB_FirePlex_chronic_irr_20190620.csv" 
dat <-  load_Fireplex_data_PhD(filename = url(url_file), threshold = 2.5)

########################
# Checking assumptions #
########################

# check for outliers 
outl <- dat %>% 
  group_by(cell_line, miRNA) %>%
  identify_outliers(log_exp) 
any(outl$is.extreme)

# normality assumption
model <- lm(log_exp~cell_line*miRNA, data = dat)
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# 
ggqqplot(dat, "log_exp", ggtheme = theme_bw()) +
  facet_grid(cell_line  ~ miRNA, labeller = "label_both")

# homogeneity of variance assumption
dat %>% levene_test(log_exp~cell_line*miRNA)

# homogeneity of variance within each miRNA
dat %>%
  group_by(miRNA) %>%
  levene_test(log_exp~cell_line) %>% arrange(p) %>%print(n="all")


#################
#   Statistics  #
#################

# expression values are log normal distributed
# to account for heterogeneity of variances use Welch ANOVA

## ================================================
## 3-way ANOVA 
# dat %>% anova_test(log_exp~miRNA*Irradiation*cell_line, white.adjust = TRUE) # white.adjust = TRUE throws error that system is singular

# 2-way ANOVA
two_way_ANOVA <- dat %>%
  group_by(miRNA) %>%
  anova_test(log_exp ~ cell_line*Irradiation, white.adjust = TRUE) %>%
  adjust_pvalue(method = "fdr") %>%
  as_tibble()

# show interaction
two_way_interaction <- two_way_ANOVA %>% 
  filter(Effect == "cell_line:Irradiation" & p.adj < 0.05) %>%
  pull(miRNA)
two_way_ANOVA %>% 
  filter(miRNA %in% two_way_interaction)%>% print(n="all")

# show main effect
two_way_main<- two_way_ANOVA %>%
  filter(Effect == "cell_line" & p.adj < 0.05) %>%
  pull(miRNA)
two_way_ANOVA %>% 
  filter(miRNA %in% two_way_main & Effect == "cell_line")

# no effect
two_way_ANOVA %>% 
  filter(!miRNA %in% two_way_main & Effect == "cell_line")


# 1-way ANOVA 
one_way_ANOVA <- dat %>%
  group_by(miRNA, Irradiation) %>%
  anova_test(log_exp ~ cell_line, white.adjust = T) %>%
  adjust_pvalue(method = "fdr") %>% 
  as_tibble() %>%
  filter(Irradiation == "control" & p.adj < 0.05) 
  




## Post-hoc Tests ................................................................................................
# post-hoc-test for the cell line main effect
pwc_cell_line <- dat %>%
  group_by(miRNA,Irradiation) %>%
  pairwise_t_test(log_exp~cell_line, pool.sd = F) %>% 
  filter(p.adj < 0.05 & Irradiation == "control")


## Differential expression..........................................................................................
# summary statistics 
dat_summary <- dat %>% 
  mutate(expression = expression +0.1) %>%
  group_by(cell_line, Irradiation,miRNA) %>% 
  summarize(geomean = geoMean(expression,na.rm=T),
            geosd = geoSD(expression,na.rm=T))

# calculate fold changes
dat_wide <-dat_summary %>% 
  select(-geosd) %>%
  filter(Irradiation == "control") %>% 
  spread(cell_line, geomean) %>%
  ungroup()

# calculate fold changes for each combination of cell lines  
cols <- c(7:4)
FC_cell_lines <- sapply(cols, y=min(cols)-1, data = dat_wide, function(data,x,y){
  tmp <- data[,x] %>% as.matrix()/data[y:(x-1)]
  colnames(tmp) <- str_c(paste(names(dat_wide[,x]),names(dat_wide[,y:(x-1)]))) %>%
    str_replace_all(" ", "_vs_")
  return(tmp)
}) %>% 
  as.data.frame() %>%
  log2() %>%
  cbind(dat_wide,.)

# filter rows for fold-changes > 1.5  
FC_cell_lines[apply(FC_cell_lines[, -c(1:7)], MARGIN = 1, function(x) any(abs(x) > log2(1.5))), ]
