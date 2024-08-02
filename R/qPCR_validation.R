library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)
library(EnvStats)
library(gtools)
library(devtools)
library(ggsci)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# access data
url_file_val <- "https://raw.githubusercontent.com/MBender1992/PhD/Marc/Data/PhD_MB_qPCR_validation_Fireplex_20201015.csv" 

# load data
dat_qPCR <- read_csv(url(url_file_val)) %>% 
  rename(Name = Probenname) %>%
  mutate(
    cell_line = str_replace_all(.$Name,"^\\d+_","") %>%
      str_replace_all("(_|\\s+)([:alpha:]+|\\d+|)\\d*[:alpha:]*","")
    ) %>%
  mutate(
    treatment = str_replace_all(.$Name,"^\\d+_[:alpha:]{3}-([:alpha:]{2}|\\d+)(_|\\s+)","") %>%
      str_replace_all("(\\d*$|_([:alpha:]+|\\d+)\\d*[:alpha:]*)","") %>%
      str_replace_all("K0", "con"),
    treatment = tolower(treatment)
    ) %>% 
  mutate(gene_name = tolower(gene_name))

# load miR30d data (calculated separately as measured with different cycler and normalized with different HK genes)
url_miR30 <- "https://raw.githubusercontent.com/MBender1992/PhD/Marc/Data/miR-30d.rds" 
dat_miR30d <- readRDS(url(url_miR30)) %>% mutate(gene_name = tolower(gene_name))


# calculate geoMean of HK genes for each sample
dat_HK <- dat_qPCR %>% 
  group_by(Name, cell_line) %>%
  filter(gene_type == "HK") %>%
  mutate(geomean_HK = geoMean(mean_ct, na.rm=T)) %>%
  distinct(geomean_HK)

# calculate dCT values
dat_dCT <- inner_join(dat_qPCR, dat_HK) %>%
  ungroup() %>%
  filter(gene_type != "HK") %>%
  select(-gene_type) %>%
  mutate(dCT = geomean_HK - mean_ct)

## calculate ddCT
ddCT <- dat_dCT %>% 
  group_by(gene_name,cell_line) %>% 
  filter(treatment == "con") %>%
  summarize(ctrl_dCT = mean(dCT, na.rm = T)) %>%
  inner_join(dat_dCT %>% filter(treatment != "con")) %>%
  mutate(ddCT = dCT - ctrl_dCT) %>%
  ungroup() %>%
  select(-c(Messung_1, Messung_2)) %>%
  add_row(dat_miR30d) %>% 
  rename(miRNA = gene_name) 



########################
# Checking assumptions #
########################


# plot qqplots of all groups
ddCT %>% 
  ggplot(aes(sample = ddCT, color = cell_line)) +
  geom_qq() +
  geom_qq_line() + 
  facet_grid(cell_line~gene_name)


# one-sided t-test vs null
ddCT %>% group_by(cell_line,miRNA) %>% 
  t_test(ddCT~ 1, mu = 0) %>% 
  adjust_pvalue(method = "fdr") %>%
  filter(p <= 0.05) %>% 
  arrange(miRNA)

ddCT %>% group_by(cell_line,miRNA) %>% 
  t_test(ddCT~ 1, mu = 0) %>% 
  adjust_pvalue(method = "fdr") %>%
  filter(p.adj <= 0.05) %>% 
  arrange(miRNA)

# statistics on ddCT, and show results in log space

# dCT[trt] = ref[trt] - goi[trt]
# dCT[ctrl] = ref[ctrl] - goi[ctrl]
# ddCT = dCT[ctrl] - dCT[trt]

# anders als 2008 beschrieben (ist intuitiver da negative dCT Werte mit weniger Expression und positive dCT mit einer 
# h?heren Expression korrelieren)






# extract fold changes for each miRNA and cell line from Fireplex analysis
logFC_fireplex <- folds %>% # from script "fireplex_differential_expression_Irradiation
  filter(miRNA %in% ddCT$miRNA) %>%
  select(-FC) %>%
  mutate(cell_line = str_replace_all(cell_line, "_","-")) %>%
  setNames(c("cell_line","miRNA", "logFC_Fireplex")) 

# join results from qPCR and Fireplex analysis
dat_plot <- ddCT %>%
  left_join(logFC_fireplex) %>%
  filter(!is.na(ddCT)) %>%
  mutate(miRNA = factor(miRNA, levels = mixedsort(unique(.$miRNA), decreasing = T)))


# calculation of sd on dCT values and calculation of ddCT by error propagation for ddCT values
svg("Results/qPCR_val.svg", width=8.5, height=7)
dat_plot %>% 
  filter(miRNA != "mir-135b-5p") %>% 
  ggplot(aes(cell_line,ddCT,fill = cell_line)) +
  stat_boxplot(
    geom ='errorbar',position = position_dodge(0.48),
    size=0.5,
    width = 0.2
  ) +
  geom_boxplot(
    color = "black", position = position_dodge(0.48), 
    mapping = aes(x = cell_line,y =ddCT),
    width=0.4
  ) +
  geom_point(aes(cell_line,logFC_Fireplex), fill = "red", shape = 23) +
  facet_wrap(~miRNA, scales = "free", nrow = 3) +
  geom_hline(yintercept = 0, lty = 3) +
  # theme_minimal() + 
  theme_PhD(axis.text.size = 12) +
  theme(legend.background = element_rect(colour = 'grey', fill = 'white', linetype='solid'),
        legend.position = "bottom",
        axis.line.y.left   = element_line(color = 'black'),
        axis.ticks.y = element_line(), 
        panel.grid  = element_blank(), 
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())  +
   scale_y_continuous(limits = c(-2.2,2.2),breaks = seq(-2,2, 1), guide = "axis_minor") +
  #scale_fill_brewer(palette = "Greys") +
  scale_fill_jco() + 
  ylab("log2-FC")
dev.off()




# confidence interval plot
# png("qPCR_val.png", units="in", width=9, height=6, res=600)
# dat_plot %>% 
#   group_by(gene_name, cell_line, logFC_Fireplex) %>%
#   summarize(mean = mean(ddCT), 
#             error = confInt (ddCT)) %>%
#   mutate(lower = mean - error, 
#          upper = mean + error) %>%
#   ggplot(aes(cell_line,mean,color = cell_line)) +
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1) +
#   geom_point(shape = 19, size = 2) + 
#   geom_point(aes(cell_line,logFC_Fireplex), fill = "orange", shape = 23) +
#   facet_wrap(~gene_name, scales = "free") +
#   geom_hline(yintercept = 0, lty = 3) +
#   theme_minimal() + 
#   theme(#legend.background = element_rect(colour = 'grey', fill = 'white', linetype='solid'),
#     legend.position = "bottom",
#     axis.line.y.left   = element_line(color = 'black'),
#     axis.ticks.y = element_line(), 
#     panel.grid  = element_blank(), 
#     axis.text.x = element_blank(),
#     axis.title.x = element_blank())  +
#   scale_y_continuous(limits = c(-2.66,2.66),breaks = seq(-2.5,2.5, 0.5)) +
#   scale_color_manual(values = cols[seq(12,20,2)], name = "cell line") +
#   scale_fill_manual(values = cols[seq(12,20,2)], name = "cell line") +
#   ylab("log2-FC") 
# dev.off()


