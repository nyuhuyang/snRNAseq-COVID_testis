library(Seurat)
library(magrittr)
library(dplyr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# Need 64GB ?
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

object = readRDS(file = "data/COVID_testis_10_SCT_20220421.rds")
meta.data = readRDS(file = "output/COVID_testis_10_20220421_metadata_v2.rds")
table(rownames(object@meta.data) == rownames(meta.data))
object@meta.data = meta.data
(step = c("resolutions")[1])

if(step == "resolutions"){# 64GB
    opts = data.frame(ident = c(rep("SCT_snn_res.0.01",7),
                                rep("SCT_snn_res.0.1",13),
                                rep("SCT_snn_res.0.2",15),
                                rep("SCT_snn_res.0.5",24),
                                rep("SCT_snn_res.0.8",28),#
                                rep("SCT_snn_res.0.9",31),
                                rep("SCT_snn_res.1",31),
                                rep("SCT_snn_res.2",45),
                                rep("SCT_snn_res.3",57),
                                rep("SCT_snn_res.4",64),
                                rep("SCT_snn_res.5",77)),
                      num = c(0:6,
                              0:12,
                              0:14,
                              0:23,
                              0:27,
                              0:30,
                              0:30,
                              0:44,
                              0:56,
                              0:63,
                              0:76)
                      )

    opt = opts[args,]
    print(opt)
    #==========================
    Idents(object) = opt$ident
    
    markers = FindMarkers_UMI(object, ident.1 = opt$num,
                              group.by = opt$ident,
                              assay = "SCT",
                              #min.pct = 0.01,
                              logfc.threshold = 0.25,
                                 only.pos = T#,
                                 #test.use = "MAST",
                                 #latent.vars = "nFeature_SCT"
                              )
    markers$cluster.number = opt$num
    markers$cluster = opt$ident
    
    num = opt$num
    if(args < 10) num = paste0("0",num)
    if(args < 100) num = paste0("0",num)
    if(args < 1000) num = paste0("0",num)
    
    arg = args
    if(args < 10) arg = paste0("0",arg)
    if(args < 100) arg = paste0("0",arg)
    if(args < 1000) arg = paste0("0",arg)
    
    write.csv(markers,paste0(path,arg,"_",opt$ident,"_",num, ".csv"))
}
