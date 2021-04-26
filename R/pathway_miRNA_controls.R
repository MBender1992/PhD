# sample(1:1000, 10)
# 
# neg con
# 
# hsa-miR-4747-5p,
# hsa-mir-592,
# hsa-miR-1-3p, 
# hsa-miR-577, 
# hsa-mir-4651,
# hsa-miR-206, 
# hsa-miR-619-5p, 
# hsa-miR-183-5p, 
# hsa-miR-7851-3p,
# hsa-miR-587 


library(tidyverse)
library(devtools)

url_file <- "https://raw.githubusercontent.com/MBender1992/PhD/Marc/Data/miRBase_miRNA_list.csv" 
dat <-  read.csv(file = url(url_file), head = TRUE, sep=";") %>% 
  as_tibble() %>% 
  mutate(ID = str_replace(ID, "mir", "miR")) %>%
  mutate(ID = str_replace(ID, "-\\d$",""))

set.seed(3)
neg_con_Cluster1 <- sample(dat$ID,10)
paste(neg_con_Cluster1,  collapse = ",")


# die Gene kann man mittels mirtargetlink ermitteln
# controls 1000 mal wiederholen und Mittelwerte bilden --> pathway Analyse in R
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5768164/ oder mirIntergrator?? 