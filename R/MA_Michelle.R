## load packages
library(vegan)
library(ggvegan)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(FactoMineR)

# load data
dat_species <- read.csv2("Data/MA_Michelle_Species.csv")
dat_spec2 <- read.csv2("Data/Begleitarten.csv")

## beschreibende Statistik
spec_sums <- dat_species %>% 
  select(-c("ALAE", "REPROD", "FOOD", "SPEC", "HABIT", "HUMID", "GK", "FLYDYM", "CODE")) %>%
  mutate(sums = rowSums(., na.rm = T)) %>% .$sums

dat_species$sums <- spec_sums
n_ges <- sum(spec_sums)

# humidity
dat_species %>% group_by(HUMID) %>% summarize(sums = sum(sums), mean = sum(sums)/n_ges)

# FOOD
dat_species %>% group_by(FOOD) %>% summarize(sums = sum(sums), mean = sum(sums)/n_ges)

# Habitat
dat_species %>% group_by(HABIT) %>% summarize(sums = sum(sums), mean = sum(sums)/n_ges)

# Flugdynamik
dat_species %>% group_by(FLYDYM) %>% summarize(sums = sum(sums), mean = sum(sums)/n_ges)

# ALAE
dat_species %>% group_by(ALAE) %>% summarize(sums = sum(sums), mean = sum(sums)/n_ges)

# Reproduction
dat_species %>% group_by(REPROD) %>% summarize(sums = sum(sums), mean = sum(sums)/n_ges)



####################################################
# Barplots für Bodenbeschaffenheit

# remove Begleitarten
dat_species <- dat_species %>% filter(!SPEC %in% dat_spec2$SPEC2)

dat_spec <- dat_species %>% 
  select(-c("ALAE", "REPROD", "FOOD", "SPEC", "HABIT", "HUMID", "GK", "FLYDYM", "CODE")) %>%
  t() %>% 
  as.data.frame()
colnames(dat_spec) <- dat_species$CODE

# calculate rowwise percentages (for each "quadrant")
dat_spec <- dat_spec %>% mutate(sums = rowSums(., na.rm = T))
dat_spec <- round((dat_spec/dat_spec$sums) *100, 2)  
dat_spec[is.na(dat_spec)] <- 0


# load soil data
dat_soil <- read.csv2("Data/MA_Michelle_Soil.csv")

# group by abundance
dat_abundant <- dat_soil %>% select(LOC, Al, P, S, K, Fe)
dat_rare     <- dat_soil %>% select(LOC, Mn, Ni, Cu, Zn, Pb)

# plot abundant elements
png("Results/abundant_elements.png", units="in", width=12, height=8, res=600)
dat_abundant %>% 
  gather(key = "element", value = "value", -LOC) %>%
  mutate(value = as.numeric(value)) %>%
  filter(LOC != "") %>%
  ggbarplot(x = "LOC", y = "value", fill = "element", position = position_dodge(), 
            add = "mean_se", add.params =  list(size = 0.7, width = 0.5), palette = "jco") +
  facet_wrap(~element, scales = "free", ncol = 2) +
  ylab("Elementgehalt in den Bodenproben (mg/kg)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), 
        legend.title = element_blank()) 
dev.off()


png("Results/rare_elements.png", units="in", width=12, height=8, res=600)
dat_rare %>% 
  gather(key = "element", value = "value", -LOC) %>%
  mutate(value = as.numeric(value)) %>%
  filter(LOC != "") %>%
  ggbarplot(x = "LOC", y = "value", fill = "element", position = position_dodge(), 
            add = "mean_se", add.params =  list(size = 0.7, width = 0.5), palette = "jco") +
  facet_wrap(~element, scales = "free", ncol = 2) +
  ylab("Elementgehalt in den Bodenproben (mg/kg)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), 
        legend.title = element_blank()) 
dev.off()


# collapse data based on location
dat_soil <- dat_soil %>% select(-c(REP, DATE, SOIL, LOC)) %>%
  mutate_if(is.character,as.numeric) %>%
  mutate(LOC = dat_soil$LOC) %>%
  group_by(LOC) %>%
  summarise_all(mean) %>%
  filter(LOC != "") %>%
  left_join(dat_soil[, c("LOC", "SOIL")], by = "LOC") %>%
  filter(!duplicated(LOC)) %>%
  mutate(SOIL = factor(SOIL)) %>%
  column_to_rownames("LOC") 



#######################################################
# Principal component analysis 

# dca to determine whether cca or pca is preferable
dca.res <- decorana(dat_spec)
dca.res 

# prepare data for PCA
dat_PCA <- dat_soil %>% 
  select(-c(V, Ga, Ge, As, Br, Rb, Sr, Y, Zr, Nb, Sn, Ba, Ce, Nd, Hf, W, Tl, Th, U, SOIL)) %>%
  mutate_if(is.numeric , replace_na, replace = 0) 
res.pca <- PCA(dat_PCA,  scale.unit = T, ncp = 5, graph = F)

# plot PCA as biplot 
png("Results/PCA_soiltype.png", units="in", width=12, height=8, res=600)
fviz_pca_biplot(res.pca, geom.ind = "point", pointshape = 21, pointsize = 2.5, 
                fill.ind = dat_soil$SOIL,
                col.ind = "black",  palette = "jco",
                alpha.var = 0.3, axes = c(1, 2),
                legend.title = list(fill = "Bodenart"),
                repel = TRUE) + 
  geom_text(aes(label = rownames(dat_soil)), nudge_x = 0.2) + 
  xlab("Erklärte Varianz PC1 (42.5%)") +
  ylab("Erklärte Varianz PC2 (18.0%)") +
  theme_bw() + 
  theme(text = element_text(size = 18), legend.position = "right")
dev.off()

# only species
dat_PCA <- dat_spec %>% select(-sums)
res.pca <- PCA(dat_PCA,  scale.unit = T, ncp = 5, graph = F)

png("Results/PCA_species.png", units="in", width=12, height=8, res=600)
fviz_pca_biplot(res.pca, geom.ind = "point", pointshape = 21, pointsize = 2.5, 
                fill.ind = dat_soil[1:18,]$SOIL,
                col.ind = "black",  palette = "jco",
                alpha.var = 0.3, axes = c(1, 2),
                legend.title = list(fill = "Bodenart"),
                repel = TRUE) + 
  xlab("Erklärte Varianz PC1 (17.6%)") +
  ylab("Erklärte Varianz PC2 (16.5%)") +
  geom_text(aes(label = rownames(dat_soil[1:18,])), nudge_x = 0.2) + 
  theme_bw() + 
  theme(text = element_text(size = 18), legend.position = "right")
dev.off()


# define data for combined PCA plot
dat_combined <- cbind(dat_soil[1:18,], dat_spec) %>% 
  select(-c(V, Ga, Ge, As, Br, Rb, Sr, Y, Zr, Nb, Sn, Ba, Ce, Nd, Hf, W, Tl, Th, U, SOIL, sums))
res.pca <- PCA(dat_combined,  scale.unit = T, ncp = 5, graph = F)

png("Results/PCA_all.png", units="in", width=12, height=8, res=600)
fviz_pca_biplot(res.pca, geom.ind = "point", pointshape = 21, pointsize = 2.5, 
                fill.ind = dat_soil[1:18,]$SOIL,
                col.ind = "black",  palette = "jco",
                alpha.var = 0.3, axes = c(1, 2),
                legend.title = list(fill = "Bodenart"),
                repel = TRUE) + 
  xlab("Erklärte Varianz PC1 (21.4%)") +
  ylab("Erklärte Varianz PC2 (13.5%)") +
  geom_text(aes(label = rownames(dat_soil[1:18,])), nudge_x = 0.2) + 
  theme_bw() + 
  theme(text = element_text(size = 18), legend.position = "right")
dev.off()

