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


ddCT <- get_ddCT(dat_qPCR, ID = "Name", group.ctrl = "time", ct.val = "Ct_Messung_290421")

ddCT %>%
  ggplot(aes(time,ddCT,fill = factor(condition))) +
  geom_bar(stat = "summary", fun = "mean",
            color = "black", position = position_dodge(0.48),
            mapping = aes(x = time,y =ddCT),
            width=0.4
  ) + geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
                    position = position_dodge(0.48), width = 0.25, size = 0.6) +
  theme_bw() +
  scale_fill_brewer()







