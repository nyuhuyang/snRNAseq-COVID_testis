# conda activate r4.1.1
library(Seurat)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
pred1 = readRDS("output/COVID_testis_10_20220421_blue_encode_singleR_pred.rds")
pred2 = readRDS("output/COVID_testis_10_20220421_KJGrive2019_singleR_pred.rds")
pred3 = readRDS("output/COVID_testis_10_20220421_Shami2020_singleR_pred.rds")
pred4 = readRDS("output/COVID_testis_10_20220421_Shami2020_PBMC_singleR_pred.rds")

meta.data = readRDS(file = "data/COVID_testis_10_20220421_metadata.rds")
meta.data %<>% cbind("blue_encode.label"=pred1$pruned.labels)
meta.data %<>% cbind("KJGrive2019.label"=pred2$pruned.labels)
meta.data %<>% cbind("Shami2020.label"=pred3$pruned.labels)
meta.data %<>% cbind("Shami2020_PBMC.label"=pred4$pruned.labels)

meta_data = read.csv("../seurat_resources/azimuth/PBMC/GSE164378/GSE164378_sc.meta.data_5P.csv",row.names =1)
meta_data = meta_data[!duplicated(meta_data$celltype.l3),]
meta.data$KJGrive2019.label %<>% plyr::mapvalues(from = meta_data$celltype.l3,
                                                 to = meta_data$celltype.l1)
meta.data$Shami2020_PBMC.label %<>% plyr::mapvalues(from = meta_data$celltype.l3,
                                                 to = meta_data$celltype.l1)

meta.data$blue_encode.label[is.na(meta.data$blue_encode.label)]= "Unknown"
meta.data$KJGrive2019.label[is.na(meta.data$KJGrive2019.label)]= "Unknown"
meta.data$Shami2020.label[is.na(meta.data$Shami2020.label)]= "Unknown"
meta.data$Shami2020_PBMC.label[is.na(meta.data$Shami2020_PBMC.label)]= "Unknown"

saveRDS(meta.data, file = "data/COVID_testis_10_20220421_metadata.rds")
