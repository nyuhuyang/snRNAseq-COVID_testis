# conda activate r4.1.1
#library(Seurat)
#library(SeuratWrappers)
#https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html#load-settings-and-packages-1
#https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html#diffusion-map-pseudotime
library(Seurat)
library(ggplot2)
library(destiny)
library(data.table)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#=============================
object = readRDS(file = "data/COVID_testis_10_SCT_20220421.rds")
meta.data = readRDS(file = "output/COVID_testis_10_20220421_metadata_v4.rds")
if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}

#14.4 Diffusion map pseudotime
#  Prepare a counts matrix with labeled rows and columns. 
#counts <- logcounts(SCE)  # access log-transformed counts matrix
#cellLabels <- SCE$cell_type
#colnames(counts) <- cellLabels
#dm <- DiffusionMap(t(counts))
"Error in dataset_extract_doublematrix(data, vars) : 
  Data needs to be matrix, data.frame, ExpressionSet, or SingleCellExperiment"
# What happens if you run the diffusion map on the PCs? Why would one do this?
pca <- object[["pca"]]@cell.embeddings
rownames(pca) <- object$cell_type
dm <- DiffusionMap(pca)
dpt <- DPT(dm, tips = which(cellLabels %in% "Spermatogonia")[1])

object$pseudotime_PC1 <- rank(pca[,"PC_1"])
object$pseudotime_harmony1 <- rank(object[["harmony"]]@cell.embeddings[,"harmony_1"])
object$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])    # rank cells by their dpt
object$pseudotime_dpt <- rank(dpt$dpt) 
library(slingshot)
SCE < as.SingleCellExperiment(object)
sce <- slingshot(SCE, reducedDim = 'PCA')  # no clusters
colors <- rainbow(50, alpha = 1)
plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)
object$slingPseudotime_1  <- sce$slingPseudotime_1

dm_cell.embeddings <- dm@eigenvectors
rownames(dm_cell.embeddings) = colnames(object)
dpt_cell.embeddings <- dpt@branch
rownames(dpt_cell.embeddings) = colnames(object)

object[["DC"]] <- CreateDimReducObject(embeddings = dm_cell.embeddings,
                                                 key = "DC_", assay = DefaultAssay(object))
#object[["DPT"]] <- CreateDimReducObject(embeddings = dpt_cell.embeddings,
#                                       key = "DPT_", assay = DefaultAssay(object))

saveRDS(object@meta.data,file = "output/COVID_testis_10_20220421_metadata_v5.rds")
saveRDS(object, file = "data/COVID_testis_10_SCT_20220421.rds")
