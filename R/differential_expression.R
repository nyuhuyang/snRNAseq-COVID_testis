invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#==============
csv_names = paste0("SCT_snn_res.",c(0.01, 0.1, 0.2, 0.5, 0.8, 0.9, 1, 2, 3,4,5))
csv_index = list.files("output/20220621",pattern = ".csv") %>% gsub("_.*","",.) %>% as.integer()
table(1:392 %in% csv_index)
csv_names = list.files("output/20220621",pattern = ".csv")
deg_list <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(paste0("output/20220621/",csv),row.names = 1)
        tmp = tmp[tmp$p_val_adj < 0.05,]
        tmp$gene = rownames(tmp)
        tmp %<>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
        tmp$resolution = sub("^.*_SCT","SCT",csv) %>% 
                sub(".csv","",.) %>% 
                sub("SCT_snn_res","SCT-snn-res",.) %>%
                sub("_.*","",.) %>%
                sub("SCT-snn-res","SCT_snn_res",.)
        tmp
})

deg = bind_rows(deg_list)
deg %<>% filter(p_val_adj < 0.05)
deg_list = split(deg, f = deg$resolution)
write.xlsx(deg_list, file = paste0(path,"COVID_DEGs.xlsx"),
           colNames = TRUE, borders = "surrounding")

#=================


DEGs <- readxl::read_excel(paste0(path,"2022-07-12- in SCT_snn_res.0.01 .xlsx"))
DEGs %<>% filter(p_val < 0.01 & avg_log2FC > 2)
table(DEGs$cluster)
openxlsx::write.xlsx(DEGs, file =  paste0(path,"2022-07-12-DEGs in SCT_snn_res.0.01.xlsx"),
                     colNames = TRUE,row.names = T,borders = "surrounding")

DEGs <- readxl::read_excel(paste0(path,"2022-07-12- in SCT_snn_res.0.5  .xlsx"))
DEGs %<>% filter(p_val < 0.01 & avg_log2FC > 2)
DEGs$cluster %<>% as.integer()
DEGs %<>% arrange(desc(avg_log2FC))
table(DEGs$cluster)
openxlsx::write.xlsx(DEGs, file =  paste0(path,"2022-07-12-DEGs in SCT_snn_res.0.05.xlsx"),
                     colNames = TRUE,row.names = T,borders = "surrounding")

