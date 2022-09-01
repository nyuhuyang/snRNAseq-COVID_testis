library(kableExtra)
library(magrittr)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

Spermatogonia <- c("GFRA1","ZBTB16","SALL4","DMRT1","ID4", "UCHL1", "L1TD1", "FMR1", "ZBTB43", "NR6A1", "BEND4", "DMRT1", "MORC1", "DAZL", "SYCP3")
Spermatocytes <- c("DAZL","KIT","CDCA8","ID4","SYCP3","TESMIN","NXT1","SHCBP1L","AURKA")
Round_Spermatids <- c("LYZL1","ACRV1","HEMGN")
Spermatids <- c("TXNDC8","TSSK6","OAZ3","PRM2")
Smooth_muscle <- c("COL1A2","ACTA2")
B_cells <- c("CD19","MS4A1")
T_cells <- c("CD3G","CD3D","CD3E","CD4","CD8A")
NK_cells <- c("GNLY","NKG7","GZMA","NCAM1")
Monocytes <- c("CD14","VCAN","FCGR3A","FCGR3B","FCN1","S100A8")
Endothelial_cells <- c("VWF","CLDN5","CDH5","PECAM1","SELE","FLT1")
Epithelial_cells <- c("KRT19","EPCAM","LGR5","BMI1")
Sertoli_markers <- c("PTGDS","WT1","MRO","SOX9","AMH","DHH")

marker = unique(c(Spermatogonia,Spermatocytes,Round_Spermatids,Spermatids,Smooth_muscle,
                  B_cells,T_cells,NK_cells,Monocytes,Endothelial_cells,Epithelial_cells,Sertoli_markers))
marker %<>% toupper()
marker %>% kable %>% kable_styling()

#==========================
meta.data = readRDS(file = "output/COVID_testis_10_20220421_metadata_v2.rds")
annotation <- readxl::read_excel("doc/annotation_covid.xlsx")
table(meta.data$SCT_snn_res.0.5 %in% annotation$resolution0.5)
unique(meta.data$SCT_snn_res.0.5[! meta.data$SCT_snn_res.0.5 %in% annotation$resolution0.5])
meta.data$cell_type <- plyr::mapvalues(meta.data$SCT_snn_res.0.5,
                                       from = annotation$resolution0.5,
                                       to = annotation$cell_type)
saveRDS(meta.data,file = "output/COVID_testis_10_20220421_metadata_v3.rds")
#==========================
df_GSE112013 <- data.table::fread("data/GSE112013_Combined_UMI_table.txt.gz")
meta.data <- pdf_text("https://static-content.springer.com/esm/art%3A10.1038%2Fs41422-018-0099-2/MediaObjects/41422_2018_99_MOESM9_ESM.pdf")
meta.data %<>% lapply(function(x) strsplit(x, split = "\\n")[[1]])
Colnames <- strsplit(meta.data[[1]][2],split = " ")[[1]]
meta.data[[1]] = meta.data[[1]][3:length(meta.data[[1]])]
meta.data %<>% lapply(function(tmp) {
    tmp %<>% lapply(function(x) t(as.data.frame(strsplit(x,split = " ")[[1]]))) %>%
        do.call("rbind",.)
    })  %>%
    do.call("rbind",.)
meta.data %<>% as.data.frame()
colnames(meta.data) = Colnames
rownames(meta.data) = meta.data$CellID

meta.data$orig.ident = gsub("-.*","",meta.data$CellID)
table(meta.data$CellID %in% colnames(df_GSE112013))
table(colnames(df_GSE112013) %in% meta.data$CellID)
