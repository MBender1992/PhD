## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(RBiomirGS)
library(ggpubr)
library(limma)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

filename <- "Data/cSCC_tissue_miRNA_Fireplex_160316.csv"
dat_tissue <- read_csv(filename) %>% 
  mutate(group = str_replace_all(group, "tumor", "Tumor"),
         group = str_replace_all(group, "exp_skin", "Exposed skin"),
         group = str_replace_all(group, "^skin", "Unexposed skin"),) %>%
  rename(Gruppe = group) %>%
  mutate(Gruppe = factor(Gruppe, levels = c("Unexposed skin", "Exposed skin", "Tumor")))

colnames(dat_tissue) <- gsub("mir", "miR", colnames(dat_tissue))

# dat_stat <- dat_tissue %>% select(-Gruppe) %>% column_to_rownames("ID") %>% t()
# dat_stat <- ifelse(is.infinite(log(dat_stat)), 0,log(dat_stat))
# 
# group_list <- factor(
#   x = c(rep("Tumor",5),  rep("Skin", 11)),
#   levels=c("Skin", "Tumor")                     
# )
# 
# design           <- model.matrix(~0+group_list) 
# colnames(design) <- c("Skin", "Tumor")  
# contrast_matrix  <-makeContrasts(Skin-Tumor, levels=group_list)
# 
# 
# # fit limma object, apply contrasts and eBayes smoothing to calculate pvals
# fit <- lmFit(dat_stat, design=design)
# fit_contrasts<-contrasts.fit(fit,contrast_matrix)
# fit_contrasts2<-eBayes(fit_contrasts)
# 
# topTable(fit_contrasts2, coef = 1, p.value = 0.05) 

dat_tissue <- dat_tissue %>% 
  mutate(ID = str_replace_all(.$ID, "tumor_1_", "cSCC-Tumor-")) %>%
  mutate(ID = str_replace_all(.$ID, "exp_skin_5_", "Exponierte Haut-")) %>%
  mutate(ID = str_replace_all(.$ID, "skin_6_", "Naive Haut-"))

dat_scaled_tissue <- dat_tissue %>% 
  column_to_rownames("ID") %>% 
  select(-Gruppe) %>%
  scale() %>%
  t()

# set ID as rownames so colorbar works properly
dat_colorbar <- dat_tissue %>% 
  select(c(ID, Gruppe)) %>% 
  column_to_rownames("ID")

# define colors for colorbar
colorbar <- HeatmapAnnotation(
  df =dat_colorbar,
  col = list(
    Gruppe = c(
      "Unexposed skin" = "steelblue1",
      "Exposed skin" = "#0073C2FF",
      "Tumor"="#EFC000FF")
  ),
  annotation_legend_param = list(
    Gruppe = list(nrow=1)
  )
)

col_fun = colorRamp2(c(-1.8, 0, 1.8), c(c("#cc00cc", "black", "#FFFF00")))

Ht_tissue <- Heatmap(dat_scaled_tissue,col= col_fun,
              top_annotation = colorbar,
              row_split = 4,
              column_title = c("A", "B"),
              border = T,
              column_km = 2,
              column_km_repeats = 100,
              clustering_method_row = "average",
              clustering_method_columns = "average",
              clustering_distance_row = "pearson",
              clustering_distance_column = "euclidean",
              rect_gp = gpar(col = "black",lty = 1, lwd = 1),
              row_names_gp = gpar(fontsize = 10),
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(
                title = "row Z-score",
                at = seq(-2,2,by=1),
                color_bar="continuous",
                title_position ="topcenter",
                legend_direction = "horizontal",
                legend_width = unit(4, "cm")
              ))

svg("Results/cSCC_tissue_Heatmap.svg", width=4, height=8)
draw(Ht_tissue, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
dev.off()

svg("Results/cSCC_tissue_miR200_family.svg", units="in", width=5, height=3, res=600)
dat_tissue %>% select(ID, Gruppe, `hsa-miR-200b-3p`, `hsa-miR-200c-3p`, `hsa-miR-429`) %>%
  gather(miRNA, exp, -c(ID, Gruppe)) %>%
  mutate(Gruppe = ifelse(Gruppe == "Tumor", "Tumor", "Skin")) %>%
  ggplot(aes(Gruppe, exp, fill = Gruppe)) +
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.1, alpha = 0.5)+
  facet_wrap(~miRNA, scales = "free", nrow = 1)+
  ggsci::scale_fill_jco(alpha = 0.7) +
  theme_PhD(axis.text.size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  ylab("miRNA expression (a.u.)")
dev.off()

############################
#   Cluster extraction     #
############################

# extract clusters out of Heatmap object
Ht_clusters <- extract_clusters(dat_scaled_tissue, Ht_tissue, sampleName = "ID", sampleClust = "sampleCluster", geneName = "miRNA", geneClust = "miRCluster")

# change type of tissue data
dat_tissue <- dat_tissue %>% 
  gather("miRNA", "expression", -c(ID, group)) %>%
  mutate(log_exp = log(expression+1))

# construct new data frame containing cluster information
dat_clusters <- left_join(dat_tissue, Ht_clusters$sampleCluster) %>%
  left_join(Ht_clusters$miRCluster) %>%
  select(-group)

# create list for each miRCluster
ls_miRCluster <- split(dat_clusters, f = dat_clusters$miRCluster)         

# compare Cluster 
cl_1A_tissue <- as.data.frame(lapply(summary_clusters(ls_miRCluster, 1, "A"), FUN=drop_attr))

# compare Cluster 
cl_4A_tissue <- as.data.frame(lapply(summary_clusters(ls_miRCluster, 4, "A"), FUN=drop_attr))


############################
#                          #
#           KEGG           #
#                          #
############################


############################
#       Cluster 1          #
############################

n <- 50
ctrl_list_KEGG<- readRDS(file = "Data/Pathway Analysis/ctrl_list_KEGG.rds") # file is available at https://github.com/MBender1992/PhD/blob/Marc/Data/Pathway%20Analysis/ctrl_list_GO_BP.rds
# calculate bias and mean targeted genes
res_ctrl_KEGG <-pathway_ctrl_summary(ctrl_list_KEGG, n=n)

# collect mRNA targets of the miRNAs upregulated in cluster 1A
rbiomirgs_mrnascan(objTitle = "cl_1A_tissue_predicted", mir = cl_1A_tissue$miRNA, sp = "hsa", 
                   queryType = "predicted", parallelComputing = TRUE, clusterType = "PSOCK")

# calculate GS by logistic regression
rbiomirgs_logistic(objTitle = "cl_1A_tissue_predicted_mirna_mrna_iwls_KEGG",mirna_DE = cl_1A_tissue, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_1A_tissue_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c2.cp.kegg.v7.2.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")

# total number of significantly enriched pathways
sum(cl_1A_tissue_predicted_mirna_mrna_iwls_KEGG_GS$adj.p.val < 0.05)

# remove enriched pathways with a random enrichment of more than 10 %
bias_KEGG <- names(res_ctrl_KEGG$bias[res_ctrl_KEGG$bias > 0.1])
cl_1A_tissue_KEGG_plot <- cl_1A_tissue_predicted_mirna_mrna_iwls_KEGG_GS %>% filter(!GS %in% bias_KEGG)


#plot results (volcano plot)
save_volcano_plot(cl_1A_tissue_KEGG_plot, n = 15, gsLabelSize = 2, sigColour = "red")
# plot distribution of enriched gene sets
save_volcano_bar_plot(cl_1A_tissue_KEGG_plot, print.all = TRUE)
# plot top enriched gene sets
save_volcano_bar_plot(cl_1A_tissue_KEGG_plot, print.all = FALSE)





############################
#       Cluster 4          #
############################

# collect mRNA targets of the miRNAs upregulated in cluster 4A
rbiomirgs_mrnascan(objTitle = "cl_4A_tissue_predicted", mir = cl_4A_tissue$miRNA, sp = "hsa", 
                   queryType = "predicted", parallelComputing = TRUE, clusterType = "PSOCK")

# calculate GS by logistic regression
rbiomirgs_logistic(objTitle = "cl_4A_tissue_predicted_mirna_mrna_iwls_KEGG",mirna_DE = cl_4A_tissue, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_4A_tissue_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c2.cp.kegg.v7.2.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")

# total number of significantly enriched pathways
sum(cl_4A_tissue_predicted_mirna_mrna_iwls_KEGG_GS$adj.p.val < 0.05)

# remove enriched pathways with a random enrichment of more than 10 %
bias_KEGG <- names(res_ctrl_KEGG$bias[res_ctrl_KEGG$bias > 0.1])
cl_4A_tissue_KEGG_plot <- cl_4A_tissue_predicted_mirna_mrna_iwls_KEGG_GS %>% filter(!GS %in% bias_KEGG)


#plot results (volcano plot)
save_volcano_plot(cl_4A_tissue_KEGG_plot, n = 15, gsLabelSize = 2, sigColour = "red")
# plot distribution of enriched gene sets
save_volcano_bar_plot(cl_4A_tissue_KEGG_plot, print.all = TRUE)
# plot top enriched gene sets
save_volcano_bar_plot(cl_4A_tissue_KEGG_plot, print.all = FALSE)




############################
#                          #
#           GOBP           #
#                          #
############################

# ############################
# #       Cluster 1          #
# ############################
# 
# 
# n <- 50
# ctrl_list_GO_BP<- readRDS(file = "Data/Pathway Analysis/ctrl_list_GO_BP.rds") # file is available at https://github.com/MBender1992/PhD/blob/Marc/Data/Pathway%20Analysis/ctrl_list_GO_BP.rds
# # calculate bias and mean targeted genes
# res_ctrl_GO_BP <-pathway_ctrl_summary(ctrl_list_GO_BP, n=n)
# 
# # calculate GS by logistic regression
# rbiomirgs_logistic(objTitle = "cl_1A_tissue_predicted_mirna_mrna_iwls_GO_BP",mirna_DE = cl_1A_tissue, 
#                    var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_1A_tissue_predicted_mrna_entrez_list, 
#                    mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c5.go.bp.v7.2.entrez.xls", optim_method = "IWLS", 
#                    p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")
# 
# # total number of significantly enriched pathways
# sum(cl_1A_tissue_predicted_mirna_mrna_iwls_GO_BP_GS$adj.p.val < 0.05)
# 
# # remove enriched pathways with a random enrichment of more than 10 %
# bias_GO_BP <- names(res_ctrl_GO_BP$bias[res_ctrl_GO_BP$bias > 0.1])
# cl_1A_tissue_GO_BP_plot <- cl_1A_tissue_predicted_mirna_mrna_iwls_GO_BP_GS %>% filter(!GS %in% bias_GO_BP)
# 
# 
# #plot results (volcano plot)
# save_volcano_plot(cl_1A_tissue_GO_BP_plot, n = 15, gsLabelSize = 2, sigColour = "red")
# # plot distribution of enriched gene sets
# save_volcano_bar_plot(cl_1A_tissue_GO_BP_plot, print.all = TRUE)
# # plot top enriched gene sets
# save_volcano_bar_plot(cl_1A_tissue_GO_BP_plot, print.all = FALSE)





############################
#       Cluster 4          #
############################

# # calculate GS by logistic regression
# rbiomirgs_logistic(objTitle = "cl_4A_tissue_predicted_mirna_mrna_iwls_GO_BP",mirna_DE = cl_4A_tissue, 
#                    var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_4A_tissue_predicted_mrna_entrez_list, 
#                    mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c5.go.bp.v7.2.entrez.xls", optim_method = "IWLS", 
#                    p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")
# 
# # total number of significantly enriched pathways
# sum(cl_4A_tissue_predicted_mirna_mrna_iwls_GO_BP_GS$adj.p.val < 0.05)
# 
# # remove enriched pathways with a random enrichment of more than 10 %
# bias_GO_BP <- names(res_ctrl_GO_BP$bias[res_ctrl_GO_BP$bias > 0.1])
# cl_4A_tissue_GO_BP_plot <- cl_4A_tissue_predicted_mirna_mrna_iwls_GO_BP_GS %>% filter(!GS %in% bias_GO_BP)
# 
# 
# #plot results (volcano plot)
# save_volcano_plot(cl_4A_tissue_GO_BP_plot, n = 15, gsLabelSize = 2, sigColour = "red")
# # plot distribution of enriched gene sets
# save_volcano_bar_plot(cl_4A_tissue_GO_BP_plot, print.all = TRUE)
# # plot top enriched gene sets
# save_volcano_bar_plot(cl_4A_tissue_GO_BP_plot, print.all = FALSE)