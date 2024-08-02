## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(devtools)
library(RBiomirGS)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# source cluster analysis
# source("R/Fireplex_cell_line_cluster_analysis.R") # if objects not loaded into workspace, load this file

## functions to use
# save volcano plot
save_volcano_plot <- function(data,...){
  png(paste("Results/",deparse(substitute(data)),"_volcano.png", sep = ""), units="in", width=5, height=4, res=600)
  p <- rbiomirgs_volcano(gsadfm = data,topgsLabel = TRUE, xLabel = "Modelkoeffizient",...)
  print(p)
  dev.off() 
}

# save bar plot within volcano plot
save_volcano_bar_plot <- function(data,print.all = TRUE,...){
  if (print.all == TRUE){
    png(paste("Results/",deparse(substitute(data)),"_volcano_bar.png", sep = ""), units="in", width=5, height=4, res=600)
    p <- rbiomirgs_bar(gsadfm = data,signif_only = F,gs.name = F,
                       n = "all",xLabel = "Gen-Set", yLabel = "Modelkoeffizient",...)
    print(p)
  } else {
    png(paste("Results/",deparse(substitute(data)),"_bar.png", sep = ""), units="in", width=5, height=4, res=600)
    p <- rbiomirgs_bar(gsadfm = data,signif_only = TRUE,gs.name = TRUE,xLabel = "Modelkoeffizient",
                       yTxtSize = 7, n = 50,...)
    print(p)
  }
  dev.off()
}


##################################
#   Pathway analysis  (GO_BP)    #
##################################


############################
#       Controls           #
############################

# calculate GS enrichment for 50 random controls 
n <- 50
# ctrl_list_GO_BP <- GS_controls(data= ctrl_stats,rep = n, miRdata = dat_miRBase$miRNA, sample_n = 10, 
#                          gs_file = "Data/Pathway Analysis/c5.go.bp.v7.2.entrez.xls")

# to facilitate downstream analysis the ctrl_list is saved 
# saveRDS(ctrl_list_GO_BP, file = "Data/Pathway Analysis/ctrl_list_GO_BP.rds")
ctrl_list_GO_BP<- readRDS(file = "Data/Pathway Analysis/ctrl_list_GO_BP.rds") # file is available at https://github.com/MBender1992/PhD/blob/Marc/Data/Pathway%20Analysis/ctrl_list_GO_BP.rds

# calculate bias and mean targeted genes
res_ctrl_GO_BP <-pathway_ctrl_summary(ctrl_list_GO_BP, n=n)



############################
#       Cluster 1          #
############################

# calculate GS by logistic regression
rbiomirgs_logistic(objTitle = "cl_1A_predicted_mirna_mrna_iwls_GO_BP",mirna_DE = cl_1A, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_1A_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL,gs_file = "Data/Pathway Analysis/c5.go.bp.v7.2.entrez.xls",optim_method = "IWLS", 
                   p.adj = "fdr",parallelComputing = FALSE,clusterType = "PSOCK")

# number of enriched GO terms
sum(cl_1A_predicted_mirna_mrna_iwls_GO_BP_GS$adj.p.val < 0.05)

# remove enriched pathways with a random enrichment of more than 10 % and only keep pathways that converged
bias_GO_BP <- names(res_ctrl_GO_BP$bias[res_ctrl_GO_BP$bias > 0.1])
cl_1A_GO_BP_plot <- cl_1A_predicted_mirna_mrna_iwls_GO_BP_GS %>% filter(!GS %in% bias_GO_BP & converged == "Y")# number of significantly enriched pathways


#plot results (volcano plot)
save_volcano_plot(cl_1A_GO_BP_plot, n = 5, gsLabelSize = 1.8, sigColour = "#CD534CFF")
# plot distribution of enriched gene sets
save_volcano_bar_plot(cl_1A_GO_BP_plot, print.all = TRUE)
# plot top enriched gene sets
save_volcano_bar_plot(cl_1A_GO_BP_plot, print.all = FALSE)


############################
#       Cluster 2          #
############################

# calculate GS by logistic regression
rbiomirgs_logistic(objTitle = "cl_2B_predicted_mirna_mrna_iwls_GO_BP",mirna_DE = cl_2B, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_2B_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL,gs_file = "Data/Pathway Analysis/c5.go.bp.v7.2.entrez.xls",optim_method = "IWLS", 
                   p.adj = "fdr",parallelComputing = FALSE,clusterType = "PSOCK")

# number of enriched GO terms
sum(cl_2B_predicted_mirna_mrna_iwls_GO_BP_GS$adj.p.val < 0.05)

# remove enriched pathways with a random enrichment of more than 10 % and only keep pathways that converged
cl_2B_GO_BP_plot <- cl_2B_predicted_mirna_mrna_iwls_GO_BP_GS %>% filter(!GS %in% bias_GO_BP & converged == "Y")# number of significantly enriched pathways

#plot results (volcano plot)
save_volcano_plot(cl_2B_GO_BP_plot, n = 5, gsLabelSize = 1.8, sigColour = "#CD534CFF")
# plot distribution of enriched gene sets
save_volcano_bar_plot(cl_2B_GO_BP_plot, print.all = TRUE)
# plot top enriched gene sets
save_volcano_bar_plot(cl_2B_GO_BP_plot, print.all = FALSE)


############################
#       Cluster 4          #
############################


# calculate GS by logistic regression
rbiomirgs_logistic(objTitle = "cl_4C_predicted_mirna_mrna_iwls_GO_BP",mirna_DE = cl_4C, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_4C_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL,gs_file = "Data/Pathway Analysis/c5.go.bp.v7.2.entrez.xls",optim_method = "IWLS", 
                   p.adj = "fdr",parallelComputing = FALSE,clusterType = "PSOCK")

# number of enriched GO terms
sum(cl_4C_predicted_mirna_mrna_iwls_GO_BP_GS$adj.p.val < 0.05)

# remove enriched pathways with a random enrichment of more than 10 % and only keep pathways that converged
cl_4C_GO_BP_plot <- cl_4C_predicted_mirna_mrna_iwls_GO_BP_GS %>% filter(!GS %in% bias_GO_BP & converged == "Y")# number of significantly enriched pathways

#plot results (volcano plot)
save_volcano_plot(cl_4C_GO_BP_plot, n = 5, gsLabelSize = 1.8, sigColour = "#CD534CFF")
# plot distribution of enriched gene sets
save_volcano_bar_plot(cl_4C_GO_BP_plot, print.all = TRUE)
# plot top enriched gene sets
save_volcano_bar_plot(cl_4C_GO_BP_plot, print.all = FALSE)
