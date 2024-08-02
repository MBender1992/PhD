## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(rstatix)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

# source("R/Fireplex_cell_line_statistics.R") # if objects not loaded into workspace, load this file
#load data
url_file <- "https://raw.githubusercontent.com/MBender1992/PhD/Marc/Data/PhD_MB_FirePlex_chronic_irr_20190620.csv" 
dat <-  load_Fireplex_data_PhD(filename = url(url_file), threshold = 2.5)

# Heatmap  transform data for Heatmap format
dat_Heatmap <- dat %>% 
  filter(Irradiation == "control" & miRNA %in% one_way_ANOVA$miRNA) %>%
  select(-c(Irradiation, log_exp)) %>%
  spread(miRNA, expression) %>%
  mutate(cell_line = str_replace_all(cell_line, "_","-")) %>%
  rename(Zelllinie = cell_line)

colnames(dat_Heatmap) <- gsub("mir", "miR", colnames(dat_Heatmap))

# scale data for Heatmap
dat_scaled <- dat_Heatmap %>% 
  column_to_rownames("ID") %>% 
  select(-Zelllinie) %>%
  scale() %>%
  t()

# set ID as rownames so colorbar works properly
dat_colorbar <- dat_Heatmap %>% select(c(ID, Zelllinie)) %>% 
  column_to_rownames("ID")

# define colors for colorbar
colorbar <- HeatmapAnnotation(
  df =dat_colorbar,
  col = list(
    Zelllinie = c(
      "Met-1" = "#0073C2FF",
      "Met-4" = "#EFC000FF",
      "SCC-12"="#868686FF" ,
      "SCC-13"="#CD534CFF",
      "SCL-II"="#7AA6DCFF")
  ),
  annotation_legend_param = list(
    Zelllinie = list(nrow=1)
  )
)

# define color scheme for Heatmap
# col_fun = colorRamp2(c(-1.8, 0, 1.8), c(c("#6D9EC1", "white", "#E46726")))
col_fun = colorRamp2(c(-1.8, 0, 1.8), c(c("#54278f", "black", "#FFFF00")))


col_fun = colorRamp2(c(-1.8, 0, 1.8), c(c("#cc00cc", "black", "#FFFF00")))

# draw Heatmap
Ht <- Heatmap(dat_scaled,col= col_fun,
  top_annotation = colorbar,
  row_split = 4,
  column_title = c("A", "B", "C"),
  border = T,
  column_km = 3,
  column_km_repeats = 100,
  clustering_method_row = "average",
  clustering_method_columns = "average",
  clustering_distance_row = "pearson",
  clustering_distance_column = "euclidean",
  rect_gp = gpar(col = "black",lty = 1, lwd = 1),
  row_names_gp = gpar(fontsize = 10),
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "row Z-score",
    at = seq(-2,2,by=1),
    color_bar="continuous",
    title_position ="topcenter",
    legend_direction = "horizontal",
    legend_width = unit(4, "cm")
  ))

svg("Results/cell_line_Heatmap.svg", width=3.5, height=9)
draw(Ht, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
dev.off()


############################
#   Cluster extraction     #
############################

# extract clusters out of Heatmap object
Ht_clusters <- extract_clusters(dat_scaled, Ht, sampleName = "ID", sampleClust = "sampleCluster", geneName = "miRNA", geneClust = "miRCluster")

# construct new data frame containing cluster information
dat_clusters <- dat %>% mutate(miRNA = str_replace_all(miRNA,"mir", "miR"))
dat_clusters <- left_join(filter(dat_clusters, Irradiation == "control"), Ht_clusters$sampleCluster) %>%
  left_join(Ht_clusters$miRCluster) %>% select(-c(Irradiation, cell_line)) %>%
  mutate(miRNA = str_replace_all(miRNA,"let", "hsa-let"),
         miRNA = str_replace_all(miRNA,"miR", "hsa-miR"))

# create list for each miRCluster
ls_miRCluster <- split(dat_clusters, f = dat_clusters$miRCluster)         

# compare Cluster 1A to 1B and 1C, drop attributes (important for the pathway analysis as a list threw an error) and generate a data.frame
cl_1A <- as.data.frame(lapply(summary_clusters(ls_miRCluster, 1, "A"), FUN=drop_attr))
# compare Cluster2B to 2A and 2C, drop attributes and generate a data.frame
cl_2B <- as.data.frame(lapply(summary_clusters(ls_miRCluster, 2, "B"), FUN=drop_attr))
# compare Cluster 4C to 4A and 4B, drop attributes and generate a data.frame
cl_4C <- as.data.frame(lapply(summary_clusters(ls_miRCluster, 4, "C"), FUN=drop_attr))


cl_A <- as.data.frame(lapply(summary_clusters(list(dat_clusters, 1), 1, "A"), FUN=drop_attr))
cl_B <- as.data.frame(lapply(summary_clusters(list(dat_clusters, 1), 1, "B"), FUN=drop_attr))
cl_C <- as.data.frame(lapply(summary_clusters(list(dat_clusters, 1), 1, "C"), FUN=drop_attr))

library(RBiomirGS)

# 
rbiomirgs_mrnascan(objTitle = "cl_A_predicted", mir = cl_A$miRNA, sp = "hsa", 
                   queryType = "predicted", parallelComputing = TRUE, clusterType = "PSOCK")
rbiomirgs_logistic(objTitle = "cl_A_predicted_mirna_mrna_iwls_KEGG",mirna_DE = cl_A, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_A_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c2.cp.kegg.v7.2.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")

# volcano plot
save_volcano_plot(cl_A_predicted_mirna_mrna_iwls_KEGG_GS, n = 15, gsLabelSize = 2, sigColour = "red")
# volcano plot
save_volcano_bar_plot(cl_A_predicted_mirna_mrna_iwls_KEGG_GS, print.all = TRUE)
# plot top enriched gene sets
save_volcano_bar_plot(cl_A_predicted_mirna_mrna_iwls_KEGG_GS, print.all = FALSE)




rbiomirgs_mrnascan(objTitle = "cl_B_predicted", mir = cl_B$miRNA, sp = "hsa", 
                   queryType = "predicted", parallelComputing = TRUE, clusterType = "PSOCK")
rbiomirgs_logistic(objTitle = "cl_B_predicted_mirna_mrna_iwls_KEGG",mirna_DE = cl_B, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_B_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c2.cp.kegg.v7.2.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")

# volcano plot
save_volcano_plot(cl_B_predicted_mirna_mrna_iwls_KEGG_GS, n = 15, gsLabelSize = 2, sigColour = "red")
# volcano plot
save_volcano_bar_plot(cl_B_predicted_mirna_mrna_iwls_KEGG_GS, print.all = TRUE)
# plot top enriched gene sets
save_volcano_bar_plot(cl_B_predicted_mirna_mrna_iwls_KEGG_GS, print.all = FALSE)





rbiomirgs_mrnascan(objTitle = "cl_C_predicted", mir = cl_C$miRNA, sp = "hsa", 
                   queryType = "predicted", parallelComputing = TRUE, clusterType = "PSOCK")
rbiomirgs_logistic(objTitle = "cl_C_predicted_mirna_mrna_iwls_KEGG",mirna_DE = cl_C, 
                   var_mirnaName = "miRNA",var_mirnaFC = "FC",var_mirnaP = "pvalue", mrnalist = cl_C_predicted_mrna_entrez_list, 
                   mrna_Weight = NULL, gs_file = "Data/Pathway Analysis/c2.cp.kegg.v7.2.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = FALSE, clusterType = "PSOCK")

# volcano plot
save_volcano_plot(cl_C_predicted_mirna_mrna_iwls_KEGG_GS, n = 15, gsLabelSize = 2, sigColour = "red")
# volcano plot
save_volcano_bar_plot(cl_C_predicted_mirna_mrna_iwls_KEGG_GS, print.all = TRUE)
# plot top enriched gene sets
save_volcano_bar_plot(cl_C_predicted_mirna_mrna_iwls_KEGG_GS, print.all = FALSE)