## <<<<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(devtools)
library(EnvStats)


# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# load and transform data
dat_qPCR <- read_csv("Data/PhD_MB_qPCR_transfection_mir_205_20210427.csv")%>% 
  rename(Name = Probenname) %>%
  mutate(time = str_extract(.$Name, "\\d+h")) %>%
  mutate(treatment = str_replace_all(.$Name, "_.+$","")) %>%
  mutate(condition = str_extract(.$Name, "\\d+nM.+microL"))

dat_mimic <- dat_qPCR %>% filter(type == "mimic")
dat_antago <- dat_qPCR %>% filter(type == "antago")

# ============================
## mimic 
ddCT_mimic <- get_ddCT(dat_mimic, ID = "Name", group.ctrl = "time", ct.val = "mean_ct")

ddCT_mimic %>%
  ggplot(aes(time,ddCT,fill = factor(condition))) +
  geom_bar(stat = "summary", fun = "mean",
            color = "black", position = position_dodge(0.48),
            mapping = aes(x = time,y =ddCT),
            width=0.4
  ) + geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
                    position = position_dodge(0.48), width = 0.25, size = 0.6) +
  theme_bw() +
  scale_fill_brewer()


# ============================
## antago mir 
ddCT_antago <- get_ddCT(dat_antago, ID = "Name", group.ctrl = "time", ct.val = "mean_ct")

ddCT_antago %>%
  ggplot(aes(time, ddCT, fill = factor(condition))) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
                  position = position_dodge(0.48), width = 0.25, size = 0.6)+
  geom_bar(stat = "summary", fun = "mean",
           color = "black", position = position_dodge(0.48),
           mapping = aes(x = time,y = ddCT),
           width=0.4
  )  + 
  geom_hline(yintercept= 0) +
  theme_PhD(axis.text.size = 8) +
  theme(
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    strip.background = element_blank(),
    legend.key.size = (unit(0.7,"line")),
    legend.position = "bottom"
  ) +
  scale_fill_brewer(palette = "Greys")



