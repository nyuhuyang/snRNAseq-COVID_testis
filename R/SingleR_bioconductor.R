#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.1.1 linux
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# ====== load single cell =============

object = readRDS("data/COVID_testis_10_SCT_20220421.rds")
sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                            colData=DataFrame(object@meta.data))
rm(object);GC()

references = c("blue_encode","KJGrive2019+PBMC","Shami2020","GuoJ2018","GuoJ2018+PBMC")[4]


# ====== load reference =============
if(references == "blue_encode"){
    blue_encode <- BlueprintEncodeData()
    rownames(blue_encode) = Hmisc::capitalize(tolower(rownames(blue_encode)))
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(blue_encode)
    ))
    length(common)
    table(blue_encode$label.fine)
    system.time(trained <- trainSingleR(ref = blue_encode[common,],
                                        labels=blue_encode$label.fine))
    system.time(pred <- classifySingleR(sce[common,], trained))
    saveRDS(object = pred, file = "output/COVID_testis_10_20220421_blue_encode_singleR_pred.rds")
}

if(references == "KJGrive2019+PBMC"){
    (load("../scRNAseq-SSC/data/SSCs_20180926.Rda"))
    SSCs = UpdateSeuratObject(SSCs)
    SSCs$celltype.l3 = SSCs$Cell.Types
    
    SSC_sce <- SingleCellExperiment(list(logcounts=SSCs[["RNA"]]@data),
                                    colData=DataFrame(SSCs[["celltype.l3"]]))
    rm(SSCs);GC()
    rownames(SSC_sce) %<>% toupper()
    
    # ======= load azimuth PBMC data ==============================
    path = "../seurat_resources/azimuth/PBMC/"
    counts <- Read10X(paste0(path, "GSE164378/GSM5008740_RNA_5P"))
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    meta.data = read.csv(paste0(path,"GSE164378/GSE164378_sc.meta.data_5P.csv"),row.names =1)
    table(rownames(meta.data) == colnames(counts))
    PBMC <- SingleCellExperiment(list(logcounts=log1p(t(t(counts)/size.factors))),
                                 colData=DataFrame(meta.data$celltype.l3))
    rm(counts);GC()
    PBMC$celltype.l3 = PBMC$meta.data.celltype.l3
    PBMC$meta.data.celltype.l3 = NULL
    # ====== conbime data =============
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(SSC_sce),
                                     rownames(PBMC)
    ))
    length(common)
    SSC_PBMC = cbind(SSC_sce[common,],PBMC[common,])

    system.time(trained <- trainSingleR(ref = SSC_PBMC[common,],
                                        labels=SSC_PBMC$celltype.l3))
    system.time(pred <- classifySingleR(sce[common,], trained))
    saveRDS(object = pred, file = "output/COVID_testis_10_20220421_KJGrive2019_singleR_pred.rds")
}


if(references == "Shami2020"){
    meta.data = data.table::fread("data/GSE142585/GSE142585_MergedHumanTestis4_PerCellAttributes.txt.gz")
    meta.data %<>% as.data.frame() %>% tibble::column_to_rownames("V1")

    counts = data.table::fread("data/GSE142585/GSE142585_MergedHumanTestis4_DGE.txt.gz")
    counts %<>% as.data.frame() %>% tibble::column_to_rownames("V1") 
    counts %<>% as.matrix %>% Matrix::Matrix(sparse = TRUE)

    table(colnames(counts) == rownames(meta.data))
    SSCs = CreateSeuratObject(counts,min.cells = 0,names.delim = "-",min.features = 0,meta.data = meta.data)
    SSCs %<>% NormalizeData()
    SSC_sce <- SingleCellExperiment(list(logcounts=SSCs[["RNA"]]@data),
                                    colData=DataFrame(SSCs@meta.data))
    rm(SSCs);GC()
    
    rownames(SSC_sce) %<>% toupper()
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(SSC_sce)
    ))
    length(common)
    
    table(SSC_sce$CellType)
    system.time(trained <- trainSingleR(ref = SSC_sce[common,],
                                        labels=SSC_sce$CellType))
    system.time(pred <- classifySingleR(sce[common,], trained))
    saveRDS(object = pred, file = "output/COVID_testis_10_20220421_Shami2020_singleR_pred.rds")
}

if(references == "Shami2020+PBMC"){
    meta.data = data.table::fread("data/GSE142585/GSE142585_MergedHumanTestis4_PerCellAttributes.txt.gz")
    meta.data %<>% as.data.frame() %>% tibble::column_to_rownames("V1")
    
    counts = data.table::fread("data/GSE142585/GSE142585_MergedHumanTestis4_DGE.txt.gz")
    counts %<>% as.data.frame() %>% tibble::column_to_rownames("V1") 
    counts %<>% as.matrix %>% Matrix::Matrix(sparse = TRUE)
    
    table(colnames(counts) == rownames(meta.data))
    SSCs = CreateSeuratObject(counts,min.cells = 0,names.delim = "-",min.features = 0,meta.data = meta.data)
    SSCs %<>% NormalizeData()
    SSCs %<>% subset(CellType %in% c("Macrophage", "Myoid", "Tcell"), invert = T)
    SSCs$celltype.l3 = SSCs$CellType
    SSC_sce <- SingleCellExperiment(list(logcounts=SSCs[["RNA"]]@data),
                                    colData=DataFrame(SSCs$celltype.l3))
    SSC_sce$celltype.l3 = SSC_sce$SSCs.celltype.l3
    SSC_sce$SSCs.celltype.l3 = NULL
    
    rm(SSCs);GC()
    
    rownames(SSC_sce) %<>% toupper()
    
    # ======= load azimuth PBMC data ==============================
    path = "../seurat_resources/azimuth/PBMC/"
    counts <- Read10X(paste0(path, "GSE164378/GSM5008740_RNA_5P"))
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    meta.data = read.csv(paste0(path,"GSE164378/GSE164378_sc.meta.data_5P.csv"),row.names =1)
    table(rownames(meta.data) == colnames(counts))
    PBMC <- SingleCellExperiment(list(logcounts=log1p(t(t(counts)/size.factors))),
                                 colData=DataFrame(meta.data$celltype.l3))
    rm(counts);GC()
    PBMC$celltype.l3 = PBMC$meta.data.celltype.l3
    PBMC$meta.data.celltype.l3 = NULL
    # ====== conbime data =============
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(SSC_sce),
                                     rownames(PBMC)
    ))
    length(common)
    SSC_PBMC = cbind(SSC_sce[common,],PBMC[common,])
    
    system.time(trained <- trainSingleR(ref = SSC_PBMC[common,],
                                        labels=SSC_PBMC$celltype.l3))
    system.time(pred <- classifySingleR(sce[common,], trained))
    saveRDS(object = pred, file = "output/COVID_testis_10_20220421_Shami2020_PBMC_singleR_pred.rds")
}

if(references == "GuoJ2018"){
    counts = data.table::fread("data/GSE112013/GSE112013_Combined_UMI_table.txt.gz")
    counts %<>% as.data.frame() %>% tibble::column_to_rownames("Gene") 
    counts %<>% as.matrix %>% Matrix::Matrix(sparse = TRUE)
    
    library(pdftools)
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
    
    meta.data$Final_clusters %<>% as.integer()
    meta.data$celltype.l2 <- plyr::mapvalues(meta.data$Final_clusters,
                                             from = 1:13,
                                             to = c("SSCs",
                                                    "Differentiating Spermatogonia",
                                                    "Early primary Spermatocyes",
                                                    "Late primary Spermatocyes",
                                                    "Round Spermatotids",
                                                    "Elongated Spermatotids",
                                                    "Sperm",
                                                    "Sperm",
                                                    "Macrophages",
                                                    "Endothelial cells",
                                                    "Myloid sells",
                                                    "Sertoli cells",
                                                    "Leydig cells"
                                                    ))
    #saveRDS(meta.data, "data/GSE112013/41422_2018_99_MOESM9_ESM_meta.data.rds")
    meta.data = readRDS(meta.data, "data/GSE112013/41422_2018_99_MOESM9_ESM_meta.data.rds")
    
    table(colnames(counts) == rownames(meta.data))
    SSCs = CreateSeuratObject(counts,min.cells = 0,names.delim = "-",min.features = 0,meta.data = meta.data)
    SSCs %<>% NormalizeData()
    SSC_sce <- SingleCellExperiment(list(logcounts=SSCs[["RNA"]]@data),
                                    colData=DataFrame(SSCs@meta.data))
    rm(SSCs);GC()
    
    #rownames(SSC_sce) %<>% toupper()
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(SSC_sce)
    ))
    length(common)
    
    table(SSC_sce$celltype.l2)
    system.time(trained <- trainSingleR(ref = SSC_sce[common,],
                                        labels=SSC_sce$celltype.l2))
    system.time(pred <- classifySingleR(sce[common,], trained))
    saveRDS(object = pred, file = "output/COVID_testis_10_20220421_GuoJ2018_singleR_pred.rds")
}
if(references == "GuoJ2018+PBMC"){
    counts = data.table::fread("data/GSE112013/GSE112013_Combined_UMI_table.txt.gz")
    counts %<>% as.data.frame() %>% tibble::column_to_rownames("Gene") 
    counts %<>% as.matrix %>% Matrix::Matrix(sparse = TRUE)
    
    library(pdftools)
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
    
    meta.data$Final_clusters %<>% as.integer()
    meta.data$celltype.l2 <- plyr::mapvalues(meta.data$Final_clusters,
                                             from = 1:13,
                                             to = c("SSCs",
                                                    "Differentiating Spermatogonia",
                                                    "Early primary Spermatocyes",
                                                    "Late primary Spermatocyes",
                                                    "Round Spermatotids",
                                                    "Elongated Spermatotids",
                                                    "Sperm",
                                                    "Sperm",
                                                    "Macrophages",
                                                    "Endothelial cells",
                                                    "Myloid sells",
                                                    "Sertoli cells",
                                                    "Leydig cells"
                                             ))
    #saveRDS(meta.data, "data/GSE112013/41422_2018_99_MOESM9_ESM_meta.data.rds")
    meta.data = readRDS(meta.data, "data/GSE112013/41422_2018_99_MOESM9_ESM_meta.data.rds")
    table(colnames(counts) == rownames(meta.data))
    SSCs = CreateSeuratObject(counts,min.cells = 0,names.delim = "-",min.features = 0,meta.data = meta.data)
    SSCs %<>% NormalizeData()
    SSC_sce <- SingleCellExperiment(list(logcounts=SSCs[["RNA"]]@data),
                                    colData=DataFrame("celltype.l2" = meta.data$celltype.l2))
    rm(SSCs);GC()
    
    # ======= load azimuth PBMC data ==============================
    path = "../seurat_resources/azimuth/PBMC/"
    counts <- Read10X(paste0(path, "GSE164378/GSM5008740_RNA_5P"))
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    meta.data = read.csv(paste0(path,"GSE164378/GSE164378_sc.meta.data_5P.csv"),row.names =1)
    table(rownames(meta.data) == colnames(counts))
    PBMC <- SingleCellExperiment(list(logcounts=log1p(t(t(counts)/size.factors))),
                                 colData=DataFrame("celltype.l2" = meta.data$celltype.l2))
    rm(counts);GC()
    PBMC$celltype.l2 = PBMC$meta.data.celltype.l2
    PBMC$meta.data.celltype.l2 = NULL
    # ====== conbime data =============
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(SSC_sce),
                                     rownames(PBMC)
    ))
    length(common)
    SSC_PBMC <- cbind(SSC_sce[common,],PBMC[common,])
    table(SSC_PBMC$celltype.l2)
    system.time(trained <- trainSingleR(ref = SSC_PBMC[common,],
                                        labels=SSC_PBMC$celltype.l2))
    system.time(pred <- classifySingleR(sce[common,], trained))
    saveRDS(object = pred, file = "output/COVID_testis_10_20220421_GuoJ2018_PBMC_singleR_pred.rds")
}

