## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(devtools)
library(rstatix)
library(RBiomirGS)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# source cluster analysis
# source("R/Fireplex_differential_expression_Irradiation.R") # if objects not loaded into workspace, load this file

dat_pathway <- diff_expr_irradiation %>%
  mutate(miRNA = stringr.tools::str_prefix(miRNA, "hsa-")) %>%
  mutate(FC = 2^log2FC) %>% 
  select(miRNA, FC, p) %>%
  setNames(c("miRNA", "FC", "pvalue")) %>%
  as.data.frame()

# collect mRNA targets of the miRNAs upregulated in cluster 1A
rbiomirgs_mrnascan(objTitle = "Fireplex_Irradiation_predicted", mir = dat_pathway$miRNA, sp = "hsa", 
                   queryType = "predicted",parallelComputing = TRUE,clusterType = "PSOCK")

# calculate GS by logistic regression
rbiomirgs_logistic(objTitle = "Fireplex_Irradiation_predicted_mirna_mrna_iwls_KEGG",mirna_DE = dat_pathway, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = Fireplex_Irradiation_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c2.cp.kegg.v7.2.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")


rbiomirgs_volcano(gsadfm = Fireplex_Irradiation_predicted_mirna_mrna_iwls_KEGG_GS,topgsLabel = TRUE, xLabel = "model coefficient",
                  n = 15, gsLabelSize = 2, sigColour = "red")




