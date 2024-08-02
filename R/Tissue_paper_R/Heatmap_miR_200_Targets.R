## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(readxl)
library(ggpubr)
library(ggh4x)

## source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

## load data
dat <- read_xlsx("Data/Tissue_paper/Targets-cSCC-Tumor-Tree-200123.xlsx")

## change dCT to -dCT values
dat[,3:10] <- apply(dat[,-c(1,2)],2, function(x){-x})

## scale data
dat_Heatmap <- dat 
dat_scaled <- dat_Heatmap %>% 
  column_to_rownames("ID") %>% 
  select(-Type) %>%
  scale() %>%
  t()

## define color palette
col_fun <- colorRamp2(c(-1.8, 0, 1.8), c(c("#cc00cc", "black", "#FFFF00")))

# set ID as rownames so colorbar works properly
dat_colorbar <- dat %>% 
  select(ID,Type) %>%  
  column_to_rownames("ID")

colorbar <- HeatmapAnnotation(
  df =dat_colorbar,
  col = list(
    Type = c(
      "Skin" = "#0073C2FF",
      "Tumor" = "#EFC000FF"#,
      #"SCC-12"="#868686FF" ,
      #"SCC-13"="#CD534CFF",
      #"SCL-II"="#7AA6DCFF")
    )
  ),
  annotation_legend_param = list(
    Type = list(nrow=1)
  )
)

# draw Heatmap
Ht <- Heatmap(dat_scaled,col= col_fun,
              top_annotation = colorbar,
              row_split = 2,
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
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(
                title = "Z-score",
                at = seq(-2,2,by=1),
                color_bar="continuous",
                title_position ="topcenter",
                legend_direction = "horizontal",
                legend_width = unit(4, "cm")
              ))

svg("Results/Tissue_paper/miR_200_Targets_Heatmap.svg",  width=3.5, height=6)
draw(Ht, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
dev.off()

dat_plot <- dat %>% gather("Gene", "neg_dcT", -c(ID, Type))

dat_plot %>%
  anova_test(neg_dcT~Gene*Type)

dat_plot %>% group_by(Gene) %>%
  games_howell_test(neg_dcT~Type)

svg("Results/Tissue_paper/miR_200_Targets_Boxplot.svg",  width=7, height=4)
dat_plot %>%
  ggplot(aes(Type, neg_dcT, fill = Type)) + 
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.1, alpha = 0.5)+
  facet_wrap(~Gene, scales = "free", nrow = 2) +
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
  scale_y_continuous(guide = "axis_minor") +
  ylab("-dCT")
dev.off()


dat_ratio <- dat_plot %>% 
  filter(Gene %in% c("CDH1", "CDH2")) %>%
  spread(Gene, neg_dcT) %>%
  mutate(ratio = -(CDH2/CDH1))

dat_ratio %>%
  t_test(ratio~Type)

svg("Results/Tissue_paper/miR_200_Targets_ratio.svg",  width=2, height=2)
dat_ratio %>% ggplot(aes(Type, ratio, fill = Type)) +
  geom_lowerrorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
                position = position_dodge(0.6), width = 0.25, size = 0.6) +
  geom_bar(stat = "summary", fun = "mean", color = "black", width = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.5)+
  theme_PhD(axis.text.size = 12) +
    theme(
      axis.title.x = element_blank(),
      strip.background = element_blank(), 
      legend.position = "none",
      legend.title = element_blank()
    ) +
    ylab("Ratio of relative N-Cadherin/E-Cadherin expression")+
    geom_hline(yintercept = 0, lty = 1) +
    scale_y_continuous(guide = "axis_minor",   limits = c(-16,0), breaks = c(seq(0,-15, -3))) +
    scale_fill_grey(start = 1, end = 0.8)
dev.off()


##*****************************************
## Correlation analysis of miR-200c-3p and ZEB1

## load data
dat <- read_csv("Data/Tissue_paper/miRNA_gene_correlation.csv") %>%
  mutate(type = factor(type, levels = c("Unexposed skin", "Exposed skin", "Local cSCC", "Metastatic cSCC")))

## function to plot correlations
correlation_plot <- function(data,x, y, group, x.anno = x_max-1, x.just = -1, ...){
 
  ## calculate pvalue for correlation
  formula <- as.formula(paste(y, "~", x))
  m.interaction <- lm(formula, data = dat)
  pval <- round(anova(m.interaction)$`Pr(>F)`[1],4)
  
  ## assign labels based on input
  x_label <- paste("Log2 of",  str_replace_all(x, "_", "-"),"relative gene expression")
  y_label <- paste("Log2",  str_replace_all(y, "_", "-"),"expression")
  
  x_min <- min(data[[x]])
  x_max <- round(max(data[[x]]))
  
  ## plot data
  dat %>%
    ggplot(aes(.data[[x]], .data[[y]])) +
    geom_smooth(method = "lm") +
    stat_regline_equation(label.x=c(x_max+x.just,-4), label.y=c(12,2), size = 3.5) +
    stat_regline_equation(aes(label=..adj.rr.label..),label.x=c(x_max+x.just,-4), label.y=c(11,2), size = 3.5) +
     annotate("text", x = x.anno, y = 10, label = paste("p = ", pval, sep = ""), size = 3.5) +
    geom_point(aes(color = .data[[group]]), size = 3) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    xlab(x_label) +
    ylab(y_label) +
    scale_y_continuous(guide = "axis_minor", limits = c(-4,12), breaks = c(seq(-4,12, 2)))+
    scale_color_manual(values = c("steelblue1", "#0073C2FF", "#EAD800ff", "#EFC000FF"))
}

## vectorize miRNA input
miR200_family <- c("miR_200b_3p", "miR_200c_3p", "miR_429") 

## apply correlation plot function
ls_ZEB1 <-lapply(miR200_family, correlation_plot, data = dat, x = "ZEB1", group = "type", x.anno = -1.89, x.just = -1.2)
ls_CDH1 <-lapply(miR200_family, correlation_plot, data = dat, x = "CDH1", group = "type", x.anno = -1.89, x.just = -1.5)
ls_CDH2 <-lapply(miR200_family, correlation_plot, data = dat, x = "CDH2", group = "type", x.anno = -5.92, x.just = -1.5)
ls_FN1 <-lapply(miR200_family, correlation_plot, data = dat, x = "FN1", group = "type", x.anno = 2.85, x.just = -1.8)
ls <- c(ls_ZEB1,ls_CDH1, ls_CDH2, ls_FN1)

svg("Results/Tissue_paper/miR200_correlation.svg", width=15, height=18)
ggarrange(plotlist = ls, common.legend = TRUE, legend = "bottom", nrow = 4, ncol = 3)
dev.off()

