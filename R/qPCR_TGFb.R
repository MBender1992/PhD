## <<<<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(devtools)
library(EnvStats)
library(readxl)
library(rstatix)
library(ggh4x)
library(ggsci)

# Messung 1 von t4 (10.11.2022) war unbrauchbar (unterschiedliche Volumina wegen Multipipette)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# load and transform data
dat_qPCR <- read_csv("Data/PhD_MB_qPCR_TGFb_20221011.csv")%>% 
  rename(Name = Probenname) %>%
  mutate(time = str_extract(.$Name, "t\\d+")) %>%
  mutate(treatment = str_replace_all(.$Name, "^.+_","")) %>%
  mutate(cell_line = str_replace_all(.$Name, "_.+",""))

dat_qPCR$treatment <- gsub("[[:digit:]]+", "", dat_qPCR$treatment)

dat_HK <- dat_qPCR %>% 
  group_by(Name) %>%
  filter(gene_type == "HK") %>%
  mutate(geomean_HK = geoMean(Ct_Messung_2, na.rm=T))  %>%
  distinct(geomean_HK) 

dat_dCT <- inner_join(dat_qPCR, dat_HK, by = c("Name")) %>%
  ungroup() %>%
  filter(gene_type != "HK") %>%
  select(-gene_type) %>%
  mutate(dCT = geomean_HK - Ct_Messung_2)


ddCT <- dat_dCT %>% 
  group_by(gene_name, cell_line, time) %>%
  filter(treatment == "ctrl") %>%
  summarize(ctrl_dCT = mean(dCT, na.rm = T)) %>%
  inner_join(dat_dCT %>% filter(treatment != "ctrl")) %>%
  mutate(ddCT = dCT - ctrl_dCT) %>%
  ungroup()

svg("Results/TGFb_Pathway_UV.svg")
ddCT %>% 
  filter(!is.na(ddCT)) %>% 
  ggplot(aes(gene_name, ddCT, fill = cell_line)) +
  # geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
  #               position = position_dodge(0.48), width = 0.25, size = 0.6) + # add only errorbars to one side of the bars
  # geom_bar(stat = "summary", fun = "mean",
  #          color = "black", position = position_dodge(0.48),
  #          mapping = aes(x = gene_name,y = ddCT),
  #          width=0.4
  # )  +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  # geom_jitter(position = position_jitterdodge(0.1), alpha = 0.7) +
  facet_wrap(cell_line~time, scales = "free") +
  geom_hline(yintercept= 0, lty = 2) +
  theme_PhD(axis.text.size = 10) +
  theme(
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),
    legend.position = "none"
  ) + 
  xlab("")+
  scale_y_continuous(limits = c(-3.5, 1.5), breaks = seq(-3,1, 1), guide = "axis_minor") +
  scale_fill_jco(alpha = 0.7)
dev.off()

ddCT %>%  filter(cell_line == "FVH2819" & gene_name != "TGFBRI") %>% group_by(gene_name, time) %>% t_test(ddCT~ 1, mu = 0) %>%
  adjust_pvalue(method = "fdr")

ddCT %>%  filter(cell_line == "SCC-12") %>% group_by(gene_name, time) %>% t_test(ddCT~ 1, mu = 0) %>%
  adjust_pvalue(method = "fdr")


nested_aov <- function(dat, x, y){
  m.full <- aov(dat[[y]]~dat[[x]])
  m.red  <- aov(dat[[y]]~0)
  anova(m.red, m.full)
}

# nested ANOVA 
dat_Fibros_t4 <- ddCT %>% filter(cell_line == "FVH2819" & gene_name != "TGFBRI" & time == "t4")
dat_Fibros_t72 <- ddCT %>% filter(cell_line == "FVH2819" & gene_name != "TGFBRI" & time == "t72")
dat_SCC12_t4 <- ddCT %>% filter(cell_line == "SCC-12"  & time == "t4")
dat_SCC12_t72 <- ddCT %>% filter(cell_line == "SCC-12"  & time == "t72")


p_Fibros_t4 <- nested_aov(dat_Fibros_t4, x = "gene_name", y = "ddCT")$`Pr(>F)`[2]
p_Fibros_t72 <- nested_aov(dat_Fibros_t72, x = "gene_name", y = "ddCT")$`Pr(>F)`[2]
p_SCC12_t4 <- nested_aov(dat_SCC12_t4, x = "gene_name", y = "ddCT")$`Pr(>F)`[2]
p_SCC12_t72 <- nested_aov(dat_SCC12_t72, x = "gene_name", y = "ddCT")$`Pr(>F)`[2]



p <- c(p_Fibros_t4, p_Fibros_t72, p_SCC12_t4,p_SCC12_t72)

p.adjust(p, method = "fdr", n = length(p))

