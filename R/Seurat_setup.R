########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.1.1
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony","magrittr"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
# read sample summary list
df_samples <- readxl::read_excel("output/20220421/20220421_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
nrow(df_samples)
#======1.2 load  Seurat =========================
object = readRDS(file = "data/COVID_testis_10_20220421.rds")

table(df_samples$sample %in% object$orig.ident)
meta.data = object@meta.data
table(df_samples$sample %in% meta.data$orig.ident)
for(i in 1:length(df_samples$sample)){
    cells <- meta.data$orig.ident %in% df_samples$sample[i]
    print(df_samples$sample[i])
    print(table(cells))
    meta.data[cells,"condition"] = as.character(df_samples$condition[i])
    meta.data[cells,"Mean.Reads.per.Cell"] = df_samples$mean.reads.per.cell[i]
    meta.data[cells,"Number.of.Reads"] = df_samples$number.of.reads[i]
    meta.data[cells,"Sequencing.Saturation"] = df_samples$sequencing.saturation[i]
}
meta.data$patient = meta.data$orig.ident
meta.data$orig.ident %<>% factor(levels = df_samples$sample)
table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data

s.genes <- cc.genes$s.genes
s.genes %<>% gsub("MLF1IP","CENPU",.)
g2m.genes <- cc.genes$g2m.genes
g2m.genes %<>% plyr::mapvalues(from = c("FAM64A", "HN1"),
                               to = c("PIMREG","JPT1"))
object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
colnames(object@meta.data)[grep("Phase",colnames(object@meta.data))]="cell.cycle.phase"


# Determine the ‘dimensionality’ of the dataset  =========
npcs = 100

DefaultAssay(object) <- "RNA"
object %<>% NormalizeData()

object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = 100, verbose = FALSE)
jpeg(paste0(path,"ElbowPlot_RNA.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = 100))
dev.off()
object <- JackStraw(object, num.replicate = 20,dims = npcs)
object <- ScoreJackStraw(object, dims = 1:npcs)

for(i in 0:9){
    a = i*10+1; b = (i+1)*10
    jpeg(paste0(path,"JackStrawPlot_",a,"_",b,".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a:b))
    dev.off()
    Progress(i, 9)
}
p.values = object[["pca"]]@jackstraw@overall.p.values
ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}
ma(p.values[,"Score"])
print(npcs <- max(which(ma(p.values[,"Score"]) <=0.05)))
npcs = 70
object[["RNA"]]@scale.data = matrix(0,0,0)

#======1.6 Performing SCTransform and integration =========================

format(object.size(object),unit = "GB")
options(future.globals.maxSize= object.size(object)*10)
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
DefaultAssay(object) <- "SCT"
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(verbose = T,npcs = npcs)

jpeg(paste0(path,"S1_ElbowPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = npcs)
dev.off()

saveRDS(object, file = "data/COVID_testis_10_20220421.rds")

#======1.8 UMAP from harmony =========================
DefaultAssay(object) = "SCT"
npcs = 70
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
#system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))

object[["harmony.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                 key = "harmonyUMAP_", assay = DefaultAssay(object))
#object[["harmony.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
#                                                 key = "harmonytSNE_", assay = DefaultAssay(object))

object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
#system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))
object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
resolutions = c( 0.01, 0.1, 0.2, 0.5,0.8,0.9, 1, 2, 3, 4, 5)
for(i in 1:length(resolutions)){
    # object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    cl = paste0("SCT_snn_res.",resolutions[i])
    print(cl)
    object@meta.data[,cl] %<>% factor(levels = as.character(sort(as.integer(levels(object@meta.data[,cl])))))
    print(table(object@meta.data[,paste0("SCT_snn_res.",resolutions[i])]))
    #Progress(i,length(resolutions))
}

saveRDS(object, file = "data/COVID_testis_10_20220421.rds")
saveRDS(object@meta.data, file = "output/COVID_testis_10_20220421_metadata_v2.rds")

#=======1.9 save SCT only =======================================
format(object.size(object),unit = "GB")

format(object.size(object@assays$RNA),unit = "GB")
object[['RNA']] <- NULL
object[["SCT"]]@counts = matrix(0,0,0)
object[["SCT"]]@scale.data = matrix(0,0,0)
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/COVID_testis_10_SCT_20220421.rds")


