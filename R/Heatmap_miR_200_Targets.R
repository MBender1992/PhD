## <<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(readxl)
library(ggpubr)
library(ggh4x)

# source R functions
source_url("https://raw.githubusercontent.com/MBender1992/base_scripts/Marc/R_functions.R")  

dat <- read_xlsx("Data/Targets-cSCC-Tumor-Tree-200123.xlsx")

dat_Heatmap <- dat %>% select(-c(ID,Type)) %>% t()
colnames(dat_Heatmap) <- dat$ID

col_fun <- colorRamp2(c(-3, 0, 3), c(c("#54278f", "black", "#FFFF00")))

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
Ht <- Heatmap(dat_Heatmap,col= col_fun,
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
                title = "Log2-FC",
                at = seq(-4,4,by=2),
                color_bar="continuous",
                title_position ="topcenter",
                legend_direction = "horizontal",
                legend_width = unit(4, "cm")
              ))

png("Results/miR_200_Targets_Heatmap.png", units="in", width=3.5, height=6, res=600)
draw(Ht, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
dev.off()

dat_plot <- dat %>% gather("Gene", "log_exp", -c(ID, Type))

dat_plot %>%
  anova_test(log_exp~Gene*Type)

dat_plot %>% group_by(Gene) %>%
  games_howell_test(log_exp~Type)

png("Results/miR_200_Targets_Boxplot.png", units="in", width=7, height=4, res=600)
dat_plot %>%
  ggplot(aes(Type, log_exp, fill = Type)) + 
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
  ylab("Log2 gene expression")
dev.off()


dat_ratio <- dat_plot %>% 
  filter(Gene %in% c("CDH1", "CDH2")) %>%
  spread(Gene, log_exp) %>%
  mutate(ratio = CDH2-CDH1)


png("Results/miR_200_Targets_ratio.png", units="in", width=2, height=2, res=600)
dat_ratio %>% ggplot(aes(Type, ratio, fill = Type)) +
  geom_uperrorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1),
                position = position_dodge(0.6), width = 0.25, size = 0.6) +
  geom_bar(stat = "summary", fun = "mean", color = "black", width = 0.6) +
  theme_PhD(axis.text.size = 12) +
    theme(
      axis.title.x = element_blank(),
      # axis.ticks.x = element_blank(),
      # axis.text.x = element_blank(),
      strip.background = element_blank(), 
      legend.position = "none",
      legend.title = element_blank()
    ) +
    ylab("N-Cadherin/ \n E-Cadherin (log2)")+
    geom_hline(yintercept = 0, lty = 1) +
    scale_y_continuous(guide = "axis_minor", expand = c(0,0),  limits = c(0,12), breaks = c(seq(-0,12, 2))) +
    scale_fill_grey(start = 0.8, end = 0.2)
dev.off()
