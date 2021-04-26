## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(devtools)
library(factoextra)
library(FactoMineR)


# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

#load data
url_file <- "https://raw.githubusercontent.com/MBender1992/PhD/Marc/Data/PhD_MB_FirePlex_chronic_irr_20190620.csv" 
dat <-  load_Fireplex_data_PhD(filename = url(url_file), threshold = 2.5)

# transform data
dat_PCA <- dat %>% 
  select(-c(expression)) %>%
  spread(miRNA, log_exp) %>%
  mutate(cell_line = str_replace_all(cell_line, "_","-"))

# res.pca <- prcomp(dat_PCA[, -c(1:3)],  scale = TRUE)
res.pca <- PCA(dat_PCA[, -c(1:3)],  scale.unit = T, ncp = 5, graph = F)

# extract results and chose which PCs to use
eig.val <- get_eigenvalue(res.pca)

# print scree plot
png("Results/cell_lines_scree_plot.png", units="in", width=5, height=4, res=600)
fviz_screeplot(res.pca) 
dev.off()

# inspect the correlation between the original variables and the principal components
fviz_pca_var(res.pca, col.var = "black")

# plot PCA by cell lines
png("Results/Fireplex_cell_line_PCA.png", units="in", width=7, height=6, res=600)
fviz_pca_ind(res.pca,  geom.ind = "point", fill.ind = dat_PCA$cell_line, col.ind = "black",
             pointshape = 21,   pointsize = dat_PCA$Irradiation, palette = "jco",
             addEllipses = T,   axes = c(2, 3), alpha.var =0.5,  mean.point = F,
             ellipse.level=0.8, legend.title = list(fill = "Cell line")) + 
  scale_size_manual(values = c(2,4))
dev.off()

# plot PCA by Irradiation
png("Results/Fireplex_irradiation_PCA.png", units="in", width=7, height=6, res=600)
fviz_pca_biplot(res.pca, geom.ind = "point", pointshape = 21,pointsize = 2.5, 
                fill.ind = dat_PCA$Irradiation, col.ind = "black",  palette = "jco",
                col.var = "cos2", alpha.var = 0.3, axes = c(1, 3),addEllipses = T, ellipse.level=0.8,
                legend.title = list(fill = "Irradiation", color = "quality of representation"),
                repel = TRUE) + 
  scale_size_manual(values = c(2,4))
dev.off()