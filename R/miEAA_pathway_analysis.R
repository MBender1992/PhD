## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(devtools)


# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

filename <- "Data/miEAA_KEGG_GO.csv"
dat_miEAA <- read_csv(filename) %>%
  mutate(`Q-value` = -log10(`Q-value`)) %>%
  mutate(Subcategory = factor(Subcategory)) %>%
  mutate(Subcategory = fct_reorder(Subcategory, `Q-value`, .desc = TRUE)) %>%
  mutate(Category = factor(Category, levels = c("KEGG (miRPathDB)", "Gene Ontology (miRWalk)")))



p <- dat_miEAA %>%
  ggplot(aes(Subcategory, `Q-value`)) + 
  geom_bar(stat = "identity", fill = "steelblue3", color = "black", width = 0.6) +
  coord_flip() +
  facet_wrap(~Category, scales = "free") + 
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("") +
  ylab("-log10(Q-value)")


png("Results/miEAA_KEGG_GO.png", units="in", width=12, height=5, res=600)
print(p)
dev.off()

