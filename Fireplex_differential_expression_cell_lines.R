##Präambel
library(tidyverse)
library(ggpubr)
library(rstatix)
# library(EnvStats)
library(ComplexHeatmap)
library(factoextra)
library(FactoMineR)
library(circlize)


if(str_detect(getwd(), "C:/")){
  setwd("C:/Marc/Arbeit/R/base_scripts")
} else{
  setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/base_scripts")
}

#loading custom functions
source("R_functions.R")


# # set working directory depending on location
# if(str_detect(getwd(), "C:/Marc")){
#   setwd("C:/Marc/Arbeit/R/PhD/Daten")
# } else{
#   setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/PhD/Daten")
# }

#load data
url_file <- "https://raw.githubusercontent.com/MBender1992/PhD/Marc/Data/200619_chronic_irr_normalized.csv" 
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


