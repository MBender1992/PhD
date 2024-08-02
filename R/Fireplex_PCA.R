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
  mutate(cell_line = str_replace_all(cell_line, "_","-")) %>%
  mutate(Irradiation = ifelse(Irradiation == "KAUVIR", "cSS", "Kontrolle")) %>%
  rename(Bestrahlung = Irradiation)

dat_PCA$Bestrahlung <- factor(dat_PCA$Bestrahlung, levels = c("Kontrolle", "cSS"))

# res.pca <- prcomp(dat_PCA[, -c(1:3)],  scale = TRUE)
res.pca <- PCA(dat_PCA[, -c(1:3)],  scale.unit = T, ncp = 5, graph = F)

# extract results and chose which PCs to use
eig.val <- get_eigenvalue(res.pca)

# print scree plot
svg("Results/cell_lines_scree_plot.svg", width=4, height=4)
fviz_screeplot(res.pca) + 
  ylab("Anteil erklärter Varianz") +
  xlab("Dimensionen") + 
  theme(text = element_text(size = 14))
dev.off()

# inspect the correlation between the original variables and the principal components
fviz_pca_var(res.pca, col.var = "black")

# plot PCA by cell lines
svg("Results/Figure1/Fireplex_cell_line_PCA_PC1_PC2.svg", width=7, height=4)
fviz_pca_ind(res.pca,  geom.ind = "point", fill.ind = dat_PCA$cell_line, col.ind = "black",
             pointshape = 21,   pointsize = dat_PCA$Bestrahlung, palette = "jco",
             addEllipses = T,   axes = c(1, 2), alpha.var =0.5,  mean.point = F,
             ellipse.level=0.8, legend.title = list(size = "Bestrahlung", fill = "Zelllinie")) + 
  scale_size_manual(values = c(2,4)) +
  theme(text = element_text(size = 18))
dev.off()


# plot PCA by Irradiation
svg("Results/Fireplex_irradiation_PCA_PC2_PC3_legend.svg",  width=7, height=6)
fviz_pca_biplot(res.pca, geom.ind = "point", pointshape = 21,pointsize = 2.5, 
                fill.ind = dat_PCA$Bestrahlung, col.ind = "black",  palette = "jco",
                col.var = "cos2", alpha.var = 0.3, axes = c(2, 3), addEllipses = T, ellipse.level=0.8,
                legend.title = list(fill = "Bestrahlung", color = "Repräsentationsqualität (cos²)"),
                repel = TRUE) + 
  scale_size_manual(values = c(2,4)) +
  theme(text = element_text(size = 18), legend.position = "none")
dev.off()


