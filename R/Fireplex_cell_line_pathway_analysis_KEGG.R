## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(devtools)
library(RBiomirGS)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# source cluster analysis
# source("R/Fireplex_cell_line_cluster_analysis.R") # if objects not loaded into workspace, load this file


############################
#       Controls           #
############################
dat_miRBase <- read.csv("Data/Pathway Analysis/mature_miRNA.csv", header = FALSE) %>%
  mutate(miRNA = str_extract(V1, "(hsa-miR|hsa-let)-\\d{1,5}([:alpha:]+-\\dp|\\s|-\\dp|[:alpha:]+)")) %>% 
  na.omit()

# data frame containing FC and pvalue of the original data used for random sampling in controls
ctrl_stats <- rbind(cl_1A, cl_2B, cl_4C) %>% 
  select(-miRNA) 

# calculate GS enrichment for 50 random controls 
n <- 50
# ctrl_list_KEGG <- GS_controls(data= ctrl_stats,rep = n, miRdata = dat_miRBase$miRNA, sample_n = 10,
#                          gs_file = "Data/Pathway Analysis/c2.cp.kegg.v7.2.entrez.gmt")

# to facilitate downstream analysis the ctrl_list is saved 
# saveRDS(ctrl_list, file = "Data/Pathway Analysis/ctrl_list.rds")
ctrl_list_KEGG <- readRDS(file = "Data/Pathway Analysis/ctrl_list_KEGG.rds") # file is available at https://github.com/MBender1992/PhD/blob/Marc/Data/Pathway%20Analysis/ctrl_list_KEGG.rds

# calculate bias and mean targeted genes
res_ctrl_KEGG <- pathway_ctrl_summary(ctrl_list_KEGG, n=n)


############################
#       Cluster 1          #
############################

# collect mRNA targets of the miRNAs upregulated in cluster 1A
rbiomirgs_mrnascan(objTitle = "cl_1A_predicted", mir = cl_1A$miRNA, sp = "hsa", 
                   queryType = "predicted",parallelComputing = TRUE, clusterType = "PSOCK")

# calculate GS by logistic regression
rbiomirgs_logistic(objTitle = "cl_1A_predicted_mirna_mrna_iwls_KEGG",mirna_DE = cl_1A, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_1A_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c2.cp.kegg.v7.2.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")

# total number of significantly enriched pathways
sum(cl_1A_predicted_mirna_mrna_iwls_KEGG_GS$adj.p.val < 0.05)

# remove enriched pathways with a random enrichment of more than 10 %
bias_KEGG <- names(res_ctrl_KEGG$bias[res_ctrl_KEGG$bias > 0.1])
cl_1A_KEGG_plot <- cl_1A_predicted_mirna_mrna_iwls_KEGG_GS %>% filter(!GS %in% bias_KEGG)


#plot results (volcano plot)
save_volcano_plot(cl_1A_KEGG_plot, n = 15, gsLabelSize = 2, sigColour = "red")
# plot distribution of enriched gene sets
save_volcano_bar_plot(cl_1A_KEGG_plot, print.all = TRUE)
# plot top enriched gene sets
save_volcano_bar_plot(cl_1A_KEGG_plot, print.all = FALSE)


############################
#       Cluster 2          #
############################


# collect mRNA targets of the miRNAs upregulated in cluster 2B
rbiomirgs_mrnascan(objTitle = "cl_2B_predicted", mir = cl_2B$miRNA, sp = "hsa", 
                   queryType = "predicted",parallelComputing = TRUE,clusterType = "PSOCK")

# calculate GS by logistic regression
rbiomirgs_logistic(objTitle = "cl_2B_predicted_mirna_mrna_iwls_KEGG",mirna_DE = cl_2B, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_2B_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c2.cp.kegg.v7.2.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")

# total number of significantly enriched pathways
sum(cl_2B_predicted_mirna_mrna_iwls_KEGG_GS$adj.p.val < 0.05)

# remove enriched pathways with a random enrichment of more than 10 %
cl_2B_KEGG_plot <- cl_2B_predicted_mirna_mrna_iwls_KEGG_GS %>% filter(!GS %in% bias_KEGG)

#plot results (volcano plot)
save_volcano_plot(cl_2B_KEGG_plot, n = 15, gsLabelSize = 2, sigColour = "red")
# plot distribution of enriched gene sets
save_volcano_bar_plot(cl_2B_KEGG_plot, print.all = TRUE)
# plot top enriched gene sets
save_volcano_bar_plot(cl_2B_KEGG_plot, print.all = FALSE)



############################
#       Cluster 4          #
############################


# collect mRNA targets of the miRNAs upregulated in cluster 4C
rbiomirgs_mrnascan(objTitle = "cl_4C_predicted", mir = cl_4C$miRNA, sp = "hsa", 
                   queryType = "predicted",parallelComputing = TRUE,clusterType = "PSOCK")

# calculate GS by logistic regression
rbiomirgs_logistic(objTitle = "cl_4C_predicted_mirna_mrna_iwls_KEGG",mirna_DE = cl_4C, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_4C_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c2.cp.kegg.v7.2.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")

# total number of significantly enriched pathways
sum(cl_4C_predicted_mirna_mrna_iwls_KEGG_GS$adj.p.val < 0.05)

# remove enriched pathways with a random enrichment of more than 10 %
cl_4C_KEGG_plot <- cl_4C_predicted_mirna_mrna_iwls_KEGG_GS %>% filter(!GS %in% bias_KEGG)

#plot results (volcano plot)
save_volcano_plot(cl_4C_KEGG_plot, n = 15, gsLabelSize = 2, sigColour = "red")
# plot distribution of enriched gene sets
save_volcano_bar_plot(cl_4C_KEGG_plot, print.all = TRUE)
# plot top enriched gene sets
save_volcano_bar_plot(cl_4C_KEGG_plot, print.all = FALSE)
