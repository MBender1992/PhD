library(tidyverse)
library(ggpubr)
library(devtools)
library(EnvStats)


# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  


# load data
# dat_qPCR <- read_csv("Data/PhD_MB_qPCR_transfection_mir_205_20210311.csv") %>%
#   rename(Name = Probenname) %>%
#   mutate(transfection = str_replace_all(.$Name,"-\\d+|(_|\\s+)([:alpha:]+|\\d+)\\d*[:alpha:]*",""),
#          transfection = tolower(transfection)) %>%
#   mutate(treatment = str_replace_all(.$Name,"(^[:alpha:]+_|_\\d+[:alpha:]+_\\d)",""),
#          treatment = tolower(treatment)) %>%
#   mutate(conc_pmol = parse_number(str_replace_all(.$Name,"^[:alpha:]+|[:alpha:]+|-\\d+|_\\d$|_",""))) %>%
#   mutate(gene_name = tolower(gene_name))
# 
# dat_sum <- dat_qPCR %>%
#   group_by(transfection, treatment, conc_pmol, gene_name) %>%
#   summarize(mean_ct = mean(Ct, na.rm=T),sd(Ct, na.rm=T)) %>%
#   mutate(gene_type = ifelse(gene_name == "mir-205-5p","target", "HK"))
# 
# #
# dat_HK <- dat_sum %>%
#   filter(gene_type == "HK") %>%
#   mutate(geomean_HK = geoMean(mean_ct, na.rm=T)) %>%
#   distinct(geomean_HK)
# 
# # calculate dCT values
# dat_dCT <- inner_join(dat_sum, dat_HK) %>%
#   ungroup() %>%
#   filter(gene_type != "HK") %>%
#   select(-gene_type) %>%
#   mutate(dCT = geomean_HK - mean_ct)
# 
# ## calculate ddCT
# ddCT <- dat_dCT %>%
#   group_by(gene_name,transfection, conc_pmol) %>%
#   filter(treatment == "ctrl") %>%
#   summarize(ctrl_dCT = mean(dCT, na.rm = T)) %>%
#   inner_join(dat_dCT %>% filter(treatment != "ctrl")) %>%
#   mutate(ddCT = dCT - ctrl_dCT) %>%
#   ungroup() %>%
#   rename(miRNA = gene_name)
# 
# ddCT %>%
#   ggplot(aes(transfection,ddCT,fill = factor(conc_pmol))) +
#   geom_bar( stat = "identity",
#     color = "black", position = position_dodge(0.48),
#     mapping = aes(x = transfection,y =ddCT),
#     width=0.4
#   ) + theme_bw() +
#   scale_fill_brewer()




dat_qPCR <- read_csv("Data/PhD_MB_qPCR_transfection_mir_205_20210427.csv") %>% 
  rename(Name = Probenname) %>%
  mutate(time = str_extract(.$Name, "\\d+h")) %>% 
  mutate(treatment = str_replace_all(.$Name, "_.+$","")) %>%
  mutate(condition = str_extract(.$Name, "\\d+nM.+microL"))

# # summarize technical replicates in excel and biological replicates in R
# dat_HK <- dat_qPCR %>% 
#   group_by(Name) %>%
#   filter(gene_type == "HK") %>%
#   mutate(geomean_HK = geoMean(Ct, na.rm=T)) %>%
#   distinct(geomean_HK) 
# 
# # calculate dCT values
# dat_dCT <- inner_join(dat_qPCR, dat_HK) %>%
#   ungroup() %>%
#   filter(gene_type != "HK") %>%
#   select(-gene_type) %>%
#   mutate(dCT = geomean_HK - Ct)
# 
# ## calculate ddCT
# ddCT <- dat_dCT %>% 
#   group_by(time) %>% 
#   filter(treatment == "ctrl") %>%
#   summarize(ctrl_dCT = mean(dCT, na.rm = T)) %>%
#   inner_join(dat_dCT %>% filter(treatment != "ctrl")) %>%
#   mutate(ddCT = dCT - ctrl_dCT) %>%
#   ungroup() %>%
#   rename(miRNA = gene_name) 

ddCT <- get_ddCT(dat_qPCR, ID = "Name", group.ctrl = "time", ct.val = "Ct")

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







