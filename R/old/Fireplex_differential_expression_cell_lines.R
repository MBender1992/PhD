## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ComplexHeatmap)
library(factoextra)
library(FactoMineR)
library(EnvStats)
library(circlize)
library(devtools)
library(RBiomirGS)


# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

#load data
url_file <- "https://raw.githubusercontent.com/MBender1992/PhD/Marc/Data/PhD_MB_FirePlex_chronic_irr_20190620.csv" 
dat <-  load_Fireplex_data_PhD(filename = url(url_file), threshold = 2.5)


###################
# Data inspection #
###################

# visualize data
dat %>% 
  ggboxplot(
    x = "cell_line",
    y = "expression",
    color = "Irradiation",
    facet.by = "miRNA",
    scales = "free", 
    palette = "jco"
  )


#########
#  PCA  #
#########

dat_PCA <- dat %>% 
  # filter(miRNA %in% pwc_cell_line$miRNA) %>%
  select(-c(expression)) %>%
  spread(miRNA, log_exp) %>%
  mutate(cell_line = str_replace_all(cell_line, "_","-"))

# res.pca <- prcomp(dat_PCA[, -c(1:3)],  scale = TRUE)
res.pca <- PCA(dat_PCA[, -c(1:3)],  scale.unit = T, ncp = 5, graph = F)

# extract results and chose which PCs to use
eig.val <- get_eigenvalue(res.pca)

# print scree plot
png("test.png", units="in", width=5, height=4, res=600)
fviz_screeplot(res.pca) 
dev.off()


# 3 PCs account for 72% var explained

# inspect the correlation between the original variables and the principal components
var$coord
fviz_pca_var(res.pca, col.var = "black")

# plot PCA
png("test.png", units="in", width=7, height=6, res=600)
fviz_pca_ind(res.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = dat_PCA$cell_line,
                col.ind = "black",
                pointshape = 21, 
                pointsize = dat_PCA$Irradiation,
                palette = "jco",
                addEllipses = T,
                axes = c(2, 3),
                # Variables
                alpha.var =0.5,
                mean.point = F,
                ellipse.level=0.8,
                legend.title = list(fill = "Cell line")
) + scale_size_manual(values = c(2,4))
dev.off()

png("test.png", units="in", width=7, height=6, res=600)
fviz_pca_biplot(res.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = dat_PCA$Irradiation,
                col.ind = "black",
                palette = "jco",
                # Color variable by groups
                col.var = "cos2",
                alpha.var = 0.3,
                axes = c(1, 3),
                addEllipses = T,
                ellipse.level=0.8,
                legend.title = list(fill = "Irradiation", color = "quality of representation"),
                repel = TRUE        # Avoid label overplotting
)+ 
  scale_size_manual(values = c(2,4))
dev.off()


########################
# Checking assumptions #
########################

# check for outliers 
outl <- dat %>% 
  group_by(cell_line, miRNA) %>%
  identify_outliers(log_exp) 
any(outl$is.extreme)

# normality assumption
model <- lm(log_exp~cell_line*miRNA, data = dat)
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# 
ggqqplot(dat, "log_exp", ggtheme = theme_bw()) +
  facet_grid(cell_line  ~ miRNA, labeller = "label_both")

# homogeneity of variance assumption
dat %>%
  levene_test(log_exp~cell_line*miRNA)

# homogeneity of variance within each miRNA
dat %>% 
  group_by(miRNA) %>%
  levene_test(log_exp~cell_line) %>% arrange(p) %>%print(n="all")


#################
#   Statistics  #
#################

# expression values are log normal distributed
# to account for heterogeneity of variances use Welch ANOVA

## 3-way ANOVA .........................................................................................................
dat %>% 
  anova_test(log_exp~miRNA*Irradiation*cell_line, white.adjust = T)

# 2-way ANOVA
two_way_ANOVA <- dat %>%
  group_by(miRNA) %>%
  anova_test(log_exp ~ cell_line*Irradiation, white.adjust = T) %>%
  adjust_pvalue(method = "fdr") 

# show interaction
two_way_interaction <- two_way_ANOVA %>%
  filter(Effect == "cell_line:Irradiation" & p.adj < 0.05) %>%
  pull(miRNA)
two_way_ANOVA %>% 
  filter(miRNA %in% two_way_interaction)%>% print(n="all")

# show main effect
two_way_main<- two_way_ANOVA %>%
  filter(Effect == "cell_line" & p.adj < 0.05) %>%
  pull(miRNA)
two_way_ANOVA %>% 
  filter(miRNA %in% two_way_main & Effect == "cell_line")

# no effect
two_way_ANOVA %>% 
  filter(!miRNA %in% two_way_main & Effect == "cell_line")


# 1-way ANOVA 
one_way_ANOVA <- dat %>%
  group_by(miRNA, Irradiation) %>%
  anova_test(log_exp ~ cell_line, white.adjust = T) %>%
  adjust_pvalue(method = "fdr") %>% 
  filter(Irradiation == "control" & p.adj < 0.05) 




## Post-hoc Tests ................................................................................................
# post-hoc-test for the cell line main effect
pwc_cell_line <- dat %>%
  group_by(miRNA,Irradiation) %>%
  pairwise_t_test(log_exp~cell_line, pool.sd = F) %>% 
  filter(p.adj < 0.05 & Irradiation == "control")


## Differential expression..........................................................................................
# summary statistics 
dat_summary <- dat %>% 
  mutate(expression = expression +0.1) %>%
  group_by(cell_line, Irradiation,miRNA) %>% 
  summarize(geomean = geoMean(expression,na.rm=T),
            geosd = geoSD(expression,na.rm=T))

# calculate fold changes
dat_wide <-dat_summary %>% 
  select(-geosd) %>%
  filter(Irradiation == "control") %>% 
  spread(cell_line, geomean) %>%
  ungroup()

# calculate fold changes for each combination of cell lines  
cols <- c(7:4)
FC_cell_lines <- sapply(cols, y=min(cols)-1, data = dat_wide, function(data,x,y){
  tmp <- data[,x] %>% as.matrix()/data[y:(x-1)]
  colnames(tmp) <- str_c(paste(names(dat_wide[,x]),names(dat_wide[,y:(x-1)]))) %>%
    str_replace_all(" ", "_vs_")
  return(tmp)
}) %>% 
  as.data.frame() %>%
  log2() %>%
  cbind(dat_wide,.)

# filter rows for fold-changes > 1.5  
FC_cell_lines[apply(FC_cell_lines[, -c(1:7)], MARGIN = 1, function(x) any(abs(x) > log2(1.5))), ]




##############
#   Heatmap  #
##############

# Heatmap  transform data for Heatmap format
dat_Heatmap <- dat %>% 
  filter(Irradiation == "control" & miRNA %in% one_way_ANOVA$miRNA) %>%
  select(-c(Irradiation, log_exp)) %>%
  spread(miRNA, expression) %>%
  mutate(cell_line = str_replace_all(cell_line, "_","-"))
  
# scale data for Heatmap
dat_scaled <- dat_Heatmap %>% 
  column_to_rownames("ID") %>% 
  select(-cell_line) %>%
  scale() %>%
  t()

# set ID as rownames so colorbar works properly
dat_colorbar <- dat_Heatmap %>% select(c(ID, cell_line)) %>% 
  column_to_rownames("ID")

# define colors for colorbar
colorbar <- HeatmapAnnotation(
  df =dat_colorbar,
  col = list(
    cell_line = c(
      "Met-1" = "#0073C2FF",
      "Met-4" = "#EFC000FF",
      "SCC-12"="#868686FF" ,
      "SCC-13"="#CD534CFF",
      "SCL-II"="#7AA6DCFF")
    ),
  annotation_legend_param = list(
    cell_line = list(nrow=1)
    )
  )



col_fun = colorRamp2(c(-1.8, 0, 1.8), c(c("#6D9EC1", "white", "#E46726")))

Ht <- Heatmap(
  dat_scaled,
  col= col_fun,
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
  rect_gp = gpar(col = "white",lty = 1, lwd = 1),
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

png("test.png", units="in", width=4, height=7, res=600)
draw(Ht, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
dev.off()

############################
#   Cluster extraction     #
############################

# extract clusters out of Heatmap object
Ht_clusters <- extract_clusters(dat_scaled, Ht, sampleName = "ID", sampleClust = "sampleCluster", geneName = "miRNA", geneClust = "miRCluster")

# construct new data frame containing cluster information
dat_clusters <- left_join(filter(dat, Irradiation == "control"), Ht_clusters$sampleCluster) %>%
  left_join(Ht_clusters$miRCluster) %>% select(-c(Irradiation, cell_line)) %>%
  mutate(miRNA = str_replace_all(miRNA,"let", "hsa-let"),
         miRNA = str_replace_all(miRNA,"mir", "hsa-miR"))

# create list for each miRCluster
ls_miRCluster <- split(dat_clusters, f = dat_clusters$miRCluster)         

# compare Cluster 1A to 1B and 1C, drop attributes (important for the pathway analysis as a list threw an error) and generate a data.frame
cl_1A <- as.data.frame(lapply(summary_clusters(ls_miRCluster, 1, "A"), FUN=drop_attr))
# compare Cluster2B to 2A and 2C, drop attributes and generate a data.frame
cl_2B <- as.data.frame(lapply(summary_clusters(ls_miRCluster, 2, "B"), FUN=drop_attr))
# compare Cluster 4C to 4A and 4B, drop attributes and generate a data.frame
cl_4C <- as.data.frame(lapply(summary_clusters(ls_miRCluster, 4, "C"), FUN=drop_attr))




##################################
#   Pathway analysis  (KEGG)     #
##################################



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
  queryType = "predicted",parallelComputing = TRUE,clusterType = "PSOCK")

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
png("volcano_cl_1A_KEGG.png", units="in", width=5, height=4, res=600)
rbiomirgs_volcano(gsadfm = cl_1A_KEGG_plot,topgsLabel = TRUE,n = 15,gsLabelSize = 2,
  sigColour = "red",plotWidth = 250,plotHeight = 220,xLabel = "model coefficient")
dev.off()

# plot distribution of enriched gene sets
png("volcano_bar_dist_cl_1A_KEGG.png", units="in", width=6, height=3, res=600)
rbiomirgs_bar(gsadfm = cl_1A_KEGG_plot,signif_only = F,gs.name = F,
  n = "all",xLabel = "gene set", yLabel = "model coefficient", plotWidth = 250, plotHeight = 220)
dev.off()

# plot top enriched gene sets
png("volcano_bar_top15cl_1A_KEGG.png", units="in", width=4, height=4, res=600)
rbiomirgs_bar(gsadfm = cl_1A_KEGG_plot,signif_only = 15,gs.name = T,xLabel = "model coefficient",
              yTxtSize = 7, n = 15, plotWidth = 250, plotHeight = 220)
dev.off()



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
png("volcano_cl_2B_KEGG.png", units="in", width=5, height=4, res=600)
rbiomirgs_volcano(gsadfm = cl_2B_KEGG_plot,topgsLabel = TRUE,n = 15,gsLabelSize = 2,
                  sigColour = "red",plotWidth = 250,plotHeight = 220,xLabel = "model coefficient")
dev.off()

# plot distribution of enriched gene sets
png("volcano_bar_dist_cl_2B_KEGG.png", units="in", width=6, height=3, res=600)
rbiomirgs_bar(gsadfm = cl_2B_KEGG_plot,signif_only = F,gs.name = F,
              n = "all",xLabel = "gene set", yLabel = "model coefficient", plotWidth = 250, plotHeight = 220)
dev.off()

# plot top enriched gene sets
png("volcano_bar_top15cl_2B_KEGG.png", units="in", width=5, height=4, res=600)
rbiomirgs_bar(gsadfm = cl_2B_KEGG_plot,signif_only = 15,gs.name = T,xLabel = "model coefficient",
              yTxtSize = 7, n ="all", plotWidth = 250, plotHeight = 220)
dev.off()



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
png("volcano_cl_4C_KEGG.png", units="in", width=5, height=4, res=600)
rbiomirgs_volcano(gsadfm = cl_4C_KEGG_plot,topgsLabel = TRUE,n = 15,gsLabelSize = 2,
                  sigColour = "red",plotWidth = 250,plotHeight = 220,xLabel = "model coefficient")
dev.off()

# plot distribution of enriched gene sets
png("volcano_bar_dist_cl_4C_KEGG.png", units="in", width=6, height=3, res=600)
rbiomirgs_bar(gsadfm = cl_4C_KEGG_plot,signif_only = F,gs.name = F,
              n = "all",xLabel = "gene set", yLabel = "model coefficient", plotWidth = 250, plotHeight = 220)
dev.off()

# plot top enriched gene sets
png("volcano_bar_top15cl_4C_KEGG.png", units="in", width=5, height=4, res=600)
rbiomirgs_bar(gsadfm = cl_4C_KEGG_plot,signif_only = 15,gs.name = T,xLabel = "model coefficient",
              yTxtSize = 7, n = "all", plotWidth = 250, plotHeight = 220)
dev.off()






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
png("volcano_cl_1A_GO_BP.png", units="in", width=5, height=4.5, res=600)
rbiomirgs_volcano(gsadfm = cl_1A_GO_BP_plot,topgsLabel = TRUE,n = 5,gsLabelSize = 1.8,
                  sigColour = "#CD534CFF",plotWidth = 250,plotHeight = 220,xLabel = "model coefficient")
dev.off()

# plot distribution of enriched gene sets
png("volcano_bar_dist_cl_1A_GO_BP.png", units="in", width=6, height=3, res=600)
rbiomirgs_bar(gsadfm = cl_1A_GO_BP_plot,signif_only = F,gs.name = F,
              n = "all",xLabel = "gene set", yLabel = "model coefficient", plotWidth = 250, plotHeight = 220)
dev.off()

# plot top enriched gene sets
png("volcano_bar_top50cl_1A_GO_BP.png", units="in", width=7, height=4.5, res=600)
rbiomirgs_bar(gsadfm = cl_1A_GO_BP_plot,signif_only = 15,gs.name = T,xLabel = "model coefficient",
              yTxtSize = 7, n = 50, plotWidth = 250, plotHeight = 220)
dev.off()


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
png("volcano_cl_2B_GO_BP.png", units="in", width=5, height=4.55, res=600)
rbiomirgs_volcano(gsadfm = cl_2B_GO_BP_plot,topgsLabel = TRUE,n = 10,gsLabelSize = 1.8,
                  sigColour = "#CD534CFF",plotWidth = 250,plotHeight = 220,xLabel = "model coefficient")
dev.off()

# plot distribution of enriched gene sets
png("volcano_bar_dist_cl_2B_GO_BP.png", units="in", width=6, height=3, res=600)
rbiomirgs_bar(gsadfm = cl_2B_GO_BP_plot,signif_only = F,gs.name = F,
              n = "all",xLabel = "gene set", yLabel = "model coefficient", plotWidth = 250, plotHeight = 220)
dev.off()

# plot top enriched gene sets
png("volcano_bar_top50cl_2B_GO_BP.png", units="in", width=8, height=4.55, res=600)
rbiomirgs_bar(gsadfm = cl_2B_GO_BP_plot,signif_only = T,gs.name = T,xLabel = "model coefficient",
              yTxtSize = 7, n = 50, plotWidth = 250, plotHeight = 220)
dev.off()




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
png("volcano_cl_4C_GO_BP.png", units="in", width=5, height=4.5, res=600)
rbiomirgs_volcano(gsadfm = cl_4C_GO_BP_plot,topgsLabel = TRUE,n = 10,gsLabelSize = 1.8,
                  sigColour = "#CD534CFF",plotWidth = 250,plotHeight = 220,xLabel = "model coefficient")
dev.off()

# plot distribution of enriched gene sets
png("volcano_bar_dist_cl_4C_GO_BP.png", units="in", width=6, height=3, res=600)
rbiomirgs_bar(gsadfm = cl_4C_GO_BP_plot,signif_only = F,gs.name = F,
              n = "all",xLabel = "gene set", yLabel = "model coefficient", plotWidth = 250, plotHeight = 220)
dev.off()

# plot top enriched gene sets
png("volcano_bar_top50cl_4C_GO_BP.png", units="in", width=8, height=4.5, res=600)
rbiomirgs_bar(gsadfm = cl_4C_GO_BP_plot,signif_only = T,gs.name = T,xLabel = "model coefficient",
              yTxtSize = 7, n = 50, plotWidth = 250, plotHeight = 220)
dev.off()


# positiver model coefficient: pathways in control group st?rker inhibiert
# negativer model coefficient: pathways in "treatment" group st?rker inhibiert

# Hypoxic naked mole-rat brains use microRNA to coordinate hypometabolic fuels and neuroprotective defenses
# Paper for Guideline to report results from BiomirGS
