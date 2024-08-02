## <<<<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(devtools)
library(EnvStats)
library(rstatix)
library(ggh4x)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# load and transform data
dat_qPCR <- read_csv("Data/PhD_MB_qPCR_transfection_mir_205_20210427.csv") %>% 
  rename(Name = Probenname) %>%
  mutate(time = str_extract(.$Name, "\\d+h")) %>%
  mutate(treatment = str_replace_all(.$Name, "_.+$","")) %>%
  mutate(condition = str_extract(.$Name, "\\d+nM.+microL"))

dat_mimic <- dat_qPCR %>% filter(type == "mimic")
dat_antago <- dat_qPCR %>% filter(type == "antago")

# ============================
## mimic 
ddCT_mimic <- get_ddCT(dat_mimic, ID = "Name", group.ctrl = "time", ct.val = "mean_ct") %>%
  mutate(condition = factor(condition, levels = c("5nM_1.5microL", "5nM_3microL", "10nM_3microL")))

p_mimic <- ddCT_mimic %>%
  ggplot(aes(time,ddCT,fill = factor(condition))) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
                position = position_dodge(0.48), width = 0.25, size = 0.6) +
  geom_bar(stat = "summary", fun = "mean",
            color = "black", position = position_dodge(0.48),
            mapping = aes(x = time,y =ddCT),
            width=0.4
  ) + 
  theme_PhD(axis.text.size = 10) +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0,6), breaks = seq(0,6,1), expand = c(0,0)) +
  scale_fill_brewer(palette = "Greys")


# ============================
## antago mir 
ddCT_antago <- get_ddCT(dat_antago, ID = "Name", group.ctrl = "time", ct.val = "mean_ct", Messung = "Messung") %>%
  mutate(condition = factor(condition, levels = c("5nM_1.5microL", "5nM_3microL", "10nM_3microL", "50nM_3microL"))) %>%
  select(-c(gene_name, type))


# add dummy data
df <- data.frame(
  time = "144h",
  condition = unique(ddCT_antago$condition)[-1],
  ddCT = 0
)



# statistics
ddCT_antago %>% group_by(time, condition) %>% t_test(ddCT~1, mu = 0, alternative = "less") %>% adjust_pvalue(method = "fdr")

ddCT_antago %>% group_by(time, condition) %>% summarize(mean = mean(ddCT, na.rm = TRUE), sd = sd(ddCT, na.rm =TRUE))

svg("Results/Antago_miR-205_qPCR_conditions.svg", width=7, height=5)
bind_rows(ddCT_antago, df) %>% filter(!is.nan(ddCT)) %>%
  mutate(time = factor(time, levels = c("24h", "72h", "144h"))) %>%
  ggplot(aes(time, ddCT, fill = factor(condition))) +
  # geom_bar(stat = "summary", fun = "mean",
  #          color = "black", position = position_dodge(0.48),
  #          mapping = aes(x = time,y = ddCT),
  #          width=0.4
  # )  +
  # geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
  #               position = position_dodge(0.48), width = 0.25, size = 0.6)+
  geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(0.5)) +
  geom_jitter(alpha = 0.4, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)) +
  geom_hline(yintercept= 0, lty = 2) +
  theme_PhD(axis.text.size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  ) + 
  scale_y_continuous(limits = c(-4,0.5), breaks = seq(-4,0.5,0.5), guide = "axis_minor") +
  scale_fill_brewer(palette = "Blues")
dev.off()



