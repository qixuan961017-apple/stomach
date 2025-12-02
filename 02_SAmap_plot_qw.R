# 加载必要的包
library(ggalluvial)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

setwd('D:/Research/Rumen/05_SAMap/')
data <- read.table("D:/Research/Rumen/05_SAMap/qw_2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(data)[1] <- "source"
colnames(data)
data_long <- pivot_longer(
  data, 
  cols = colnames(data)[2:ncol(data)],  # 明确选择所有目标节点列
  names_to = "target",  # 目标节点列名
  values_to = "value"   # 流动值
)
data_long <- data_long[
  (grepl("^de", data_long$source) & grepl("^hu", data_long$target)) |  # hu -> de
    (grepl("^hu", data_long$source) & grepl("^SN", data_long$target)),   # de -> SN
]

# 提取 source 列中的信息
data_long <- data_long %>%
  mutate(
    Species_source = str_extract(source, "^de|^hu|^SN"),  # 提取 'de_deer', 'hu', 'SN_SN'
    Celltype_source = str_extract(source, "(?<=_)[A-Za-z ]+(?=\\s|$)"),  # 提取细胞类型 (如 Fibroblast)
    Cellcluster_source = str_extract(source, "(lw|zw)")  # 提取 'lw' 或 'zw'
  )

# 提取 target 列中的信息
data_long <- data_long %>%
  mutate(
    Species_target = str_extract(target, "^de|^hu|^SN"),  # 提取 'de_deer', 'hu', 'SN_SN'
    Celltype_target = str_extract(target, "(?<=_)[A-Za-z ]+(?=\\s|$)"),  # 提取细胞类型 (如 Fibroblast)
    Cellcluster_target = str_extract(target, "(lw|zw|LW|ZW)") ,  # 提取 'lw' 或 'zw'
    Cellcluster_target = str_replace_all(Cellcluster_target, c("LW" = "lw", "ZW" = "zw"))  # 替换 LW 和 ZW
  )

##修改细胞类型的名称
data_long$Celltype_source <- recode(data_long$Celltype_source,
                                    "Monocytes" = "Mon",
                                    "Monocyte" = "Mon",
                                    "Smooth Muscle cell" = "SMC",
                                    "T cell" = "T",
                                    "B cell" = "B",
                                    "Endothelial cell" = "End",
                                    "Endothelial" = "End",
                                    "Enteroendocrine cell" = "Ent",
                                    "Epithelial cell" = "Epi",
                                    "Epithelial" = "Epi",
                                    "Fibroblast" = "Fib",
                                    "Mast" = "MC",
                                    "Mast cell" = "MC")
data_long$Celltype_target <- recode(data_long$Celltype_target,
                                    "Monocytes" = "Mon",
                                    "Monocyte" = "Mon",
                                    "Smooth Muscle cell" = "SMC",
                                    "T cell" = "T",
                                    "B cell" = "B",
                                    "Endothelial cell" = "End",
                                    "Endothelial" = "End",
                                    "Enteroendocrine cell" = "Ent",
                                    "Epithelial cell" = "Epi",
                                    "Epithelial" = "Epi",
                                    "Fibroblast" = "Fib",
                                    "Mast" = "MC",
                                    "Mast cell" = "MC")

#install.packages("networkD3")
library(networkD3)
# 创建 Sankey 图
data_long$source <- paste(data_long$Celltype_source,data_long$Species_source,data_long$Cellcluster_source,sep = "|")
data_long$target <- paste(data_long$Celltype_target,data_long$Species_target,data_long$Cellcluster_target,sep = "|")

# 创建节点列表 
nodes <- data.frame(name = unique(c(as.character(data_long$source), as.character(data_long$target))))
# 创建链接列表
links <- data.frame(
  source = match(data_long$source, nodes$name) - 1,  # 转换为从 0 开始的索引
  target = match(data_long$target, nodes$name) - 1,  # 转换为从 0 开始的索引
  value = data_long$value
)
links <- links[links$value >= 0.1, ] #0.01
##color setting---------start---------------------------------------------------------
nodes<-data.frame(do.call(rbind,strsplit(nodes$name,'[|]')))
colnames(nodes)=c('name','species','stomach')
nodes$name <- paste(nodes$name, nodes$stomach, sep = '|')
nodes$group<-as.factor(nodes$name)

color_regions = 'd3.scaleOrdinal()
                .domain(["Mon|lw", "SMC|lw", "T|lw", "B|lw", "Ent|lw", "End|lw", "Epi|lw", "Fib|lw", "MC|lw", "Mon|zw", "SMC|zw", "T|zw", "B|zw", "Ent|zw", "End|zw", "Epi|zw", "Fib|zw", "MC|zw", "Mon|NA", "SMC|NA", "T|NA", "B|NA", "Ent|NA", "End|NA", "Epi|NA", "Fib|NA", "MC|NA"])
                .range(["#907FD3", "#91ae6c", "#95D0E2", "#5F7DC7", "#f77904", "#f3cb45", "#f12929", "#b8ce83", "#C7D1EE"])'


sankey<-sankeyNetwork(Links = links,
                   Nodes = nodes, 
                   Source = "source",
                   Target = "target", 
                   Value = "value", 
                   NodeID = "name",
                   NodeGroup = "group",
                   colourScale = JS(color_regions),
                   # LinkGroup = 'energy_type',
                   units = "votes",
                   height =800,
                   width =1000,
                   nodePadding = 20,
                   sinksRight =  FALSE, 
                   fontSize = 20, 
                   nodeWidth = 20
)
sankey
saveNetwork(sankey, file = "sankey_plot.html")
#install.packages("webshot")
#webshot::install_phantomjs()  # 必须安装phantomjs才能截图
webshot::webshot("sankey_plot.html", "sankey_plot.pdf")















###########作图二###################
library(ggalluvial)
library(ggplot2)
library(tidyr)
setwd('D:/Research/Rumen/05_SAMap/')
data <- read.table("D:/Research/Rumen/05_SAMap/qw_2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(data)[1] <- "source"
data_long <- pivot_longer(
  data, 
  cols = colnames(data)[2:ncol(data)],  # 明确选择所有目标节点列
  names_to = "target",  # 目标节点列名
  values_to = "value"   # 流动值
)
# data_long <- data_long[
#   grepl("^hu", data_long$source) & grepl("^de", data_long$target),   # hu -> de
# ]
# 过滤数据：根据你的要求，保留符合条件的行
data_long <- data_long[
  (grepl("^de", data_long$source) & grepl("^hu", data_long$target)) |  # hu -> de
    (grepl("^hu", data_long$source) & grepl("^SN", data_long$target)),   # de -> SN
]
# 定义颜色
sankey_colors <- c("#0072B2", "#999999", "#D55E00", "#E69F00", "#009E73", "#56B4E9", "#F0E442", "#111111", 
                   "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
                   "#66C2A5", "#FC8D62", "#8DA0CB", "#E6F598", "#FF0033", "#9C27B0", "#42A5F5", "#81C784", 
                   "#FF9800", "#C2185B", "#FFEB3B", "#8BC34A")

# 创建 ggalluvial Sankey 图
p.sankey <- ggplot(data = data_long,
                   aes(axis1 = source, axis2 = target, y = value)) +
  scale_x_discrete(limits = c("Source", "Target"), expand = c(0.15, 0.05)) +  # 定义轴
  geom_alluvium(aes(fill = source), alpha = 1, curve_type = "arctangent") +  # 绘制流动
  scale_color_manual(values = sankey_colors) +  # 设置颜色
  scale_fill_manual(values = sankey_colors) +  # 设置填充颜色
  geom_stratum(alpha = 0, color = adjustcolor("white", alpha.f = 1), size = 1.2, fill = 'white') +  # 设置层次样式
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), cex = 6) +  # 添加标签
  theme_void() +  # 去掉背景
  theme(
    legend.position = "none",  # 不显示图例
    axis.text = element_text(size = 24),  # 设置字体大小
    axis.title = element_blank(),  # 不显示轴标题
    axis.text.y = element_blank(),  # 不显示 y 轴标签
    panel.grid = element_blank()  # 去除网格线
  )
p.sankey
ggsave(p.sankey, filename = "sankey_plot.pdf", width = 20, height = 6)

