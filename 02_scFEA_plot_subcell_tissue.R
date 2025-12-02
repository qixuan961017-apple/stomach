#=========================scFEA downstream analysis===========================
#-----------------------------------------------------------------------------
#1、热图展示celltype每组中的代谢flux差异

library(Seurat)
library(dplyr)
setwd("D:/Research/Rumen/14_scFEA/")

adj_flux <- read.csv("./output/adj_flux.csv",row.names = 1)
rownames(adj_flux) <- gsub("\\.", "-", rownames(adj_flux))

adj_scRNA <- readRDS("D:/Research/Rumen/14_scFEA/merged_qw_zw.rds")
adj_scRNA$sample <- gsub("sheep[23]", "sheep", adj_scRNA$orig.ident)

tissue_labels <- sapply(adj_scRNA$orig.ident, function(x) {
  if (grepl("lw", x, ignore.case = TRUE)) return("lw") #大小写不敏感
  if (grepl("ww", x, ignore.case = TRUE)) return("ww")
  if (grepl("bw", x, ignore.case = TRUE)) return("bw")
  if (grepl("zw", x, ignore.case = TRUE)) return("zw")
})
adj_scRNA$tissue <- tissue_labels

#cell注释
adj_scRNA$group_cells <- paste0(adj_scRNA$subcell,"_",adj_scRNA$tissue)
cell_anno <- data.frame(cellid=rownames(adj_scRNA@meta.data),
                        group_cells=adj_scRNA$group_cells)

cell_anno <- cell_anno[order(cell_anno$group_cells),]
adj_flux <- adj_flux[cell_anno$cellid,]


#计算每组的代谢通量平均值
df_averages <- adj_flux %>%
  group_by(group = cell_anno$group_cells) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  select(-group)

rownames(df_averages) <- unique(cell_anno$group_cells)
df_averages <- t(df_averages)%>% as.data.frame()


df_flux <- df_averages

#删掉标准差为0的行
df_flux = df_flux[apply(df_flux, 1, function(x) sd(x)!=0),]

library(pheatmap)
df_flux[is.na(df_flux)] <- 0

pheatmap::pheatmap(df_flux, cluster_cols = F, cluster_rows = T,
                   show_rownames = F, scale = "row",
                   colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100),
                   border=F,heatmap_legend_param = list(title="Flux"),
                   fontsize_col=8,
                   treeheight_row=20)


#做一下注释----------------------------------------------------------------------
#列注释
library(stringr)
# samples <- gsub("^([^_]+_[^_]+).*", "\\1", colnames(df_flux))
samples <- gsub(".*_", "", colnames(df_flux)) 
celltypes <- sub("_.*", "", colnames(df_flux))
annotation_col <- data.frame(
  celltype = celltypes,
  sample = samples
)
row.names(annotation_col) <- colnames(df_flux)

#Module注释
human_moduleInfo <- read.csv("scFEA.human.moduleinfo.csv", header = T, row.names = 1)

#行注释
annotation_row = human_moduleInfo[rownames(df_flux),]
annotation_row  = as.data.frame(annotation_row[,c("SM_anno")])
rownames(annotation_row) = rownames(df_flux)
colnames(annotation_row) = c("SM_anno")


#注释color
cellcolor <- c('#e12e2b', '#eb7617', '#e9c449', '#b1c780', '#8ea86a','#91c8d9','#5b76b9', '#877bb7','#c4cde8')
names(cellcolor) <- c("Epithelial","Enteroendocrine","Endothelial","Fibroblast","SMC","T","B","Monocytes","Mast")

groupcolor <- c('#89cef1', '#f59295', '#e59230', '#81bf6c')
names(groupcolor) <- c("lw", "ww",  "bw", "zw")

modulecolor <- c("#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#FFFF00",
                 "#808000","#FF00FF","#FA8072","#800080","#87CEEB","#40E0D0","#5F9EA0",
                 "#008B8B","#FFE4B5","#228B22","#4682B4","#32CD32","#F0E68C","#FFFFE0",
                 "#FF6347")
names(modulecolor) <- unique(annotation_row$SM_anno)

ann_colors <- list(celltype=cellcolor, sample= groupcolor, SM_anno=modulecolor) #颜色设置

#宽数据转换成长数据
library(reshape2)
df_flux_with_genes <- df_flux %>% tibble::rownames_to_column("module")
df_long <- melt(df_flux_with_genes, id.vars = "module", variable.name = "sample", value.name = "value")
top10 <- df_long %>% group_by(sample) %>% top_n(n = 50, wt = value)


#标注感兴趣的module----------------------------------------------------------------------
#我们标注module和代谢产物
colnames(human_moduleInfo)
human_moduleInfo$module_name <- paste0(human_moduleInfo$Module_id,": ",
                                       human_moduleInfo$Compound_IN_name,
                                       "_", human_moduleInfo$Compound_OUT_name)
# write.csv(human_moduleInfo, file = 'scFEA_human_moduleInfo.csv')

df_flux <- df_flux[rownames(df_flux) %in% unique(top10$module),]
select_moduleInfo = human_moduleInfo[rownames(df_flux),]

#替换flux的行名
df_flux_new <- df_flux
annotation_row = human_moduleInfo[rownames(df_flux_new),]
rownames(df_flux_new) <- select_moduleInfo$module_name

#注释的名字也要变
row.names(annotation_col) <- colnames(df_flux_new)
rownames(annotation_row) = rownames(df_flux_new)
annotation_row <- annotation_row %>% select(SM_anno)


# Ensure celltype order in annotation_row matches cellcolor order
annotation_col$sample <- factor(annotation_col$sample, levels = names(groupcolor))
annotation_col$celltype  <- factor(annotation_col$celltype, levels = names(cellcolor))
# Verify sorting order of annotation_col
annotation_col <- annotation_col[order(annotation_col$sample), ]
annotation_col <- annotation_col[order(annotation_col$celltype), ]

df_flux_new <- df_flux_new[, order(match(colnames(df_flux_new), rownames(annotation_col)))]
# pdf("./3-anno-celltypeFlux.pdf", width = 10, height = 10)
ht1 = pheatmap::pheatmap(df_flux_new, scale = "row",show_rownames = T,#不显示行名
                         show_colnames = T,#不显示列名
                         cluster_cols = F, cluster_rows = T,
                         col = colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100),
                         annotation_col = annotation_col, #列注释信息
                         annotation_row = annotation_row,#行注释信息
                         annotation_names_row = F,#不显示行注释信息
                         annotation_names_col = F ,#不显示列注释信息
                         column_title = NULL,#不显示列标题
                         row_title = NULL,
                         fontsize_col=8,
                         treeheight_row=20,
                         heatmap_legend_param = list(title="Flux"),
                         annotation_colors = ann_colors)