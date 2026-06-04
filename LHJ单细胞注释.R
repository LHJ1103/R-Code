library(hdf5r)
library(Seurat) 
library(dplyr)
library(patchwork)
library(tidyverse)
library(reticulate)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readr)
library(stringr)

DefaultAssay(UCCRCHarmony)<-"RNA"
UCCRCHarmony <- ScaleData(UCCRCHarmony)
plot1=DimPlot(UCCRC,reduction = "umap",label = T)
plot2=DimPlot(UCCRC,reduction = "umap",group.by = "library")
table(UCCRC@meta.data$library)

markers<-FindAllMarkers(object = UCCRC,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
top30<-markers%>%group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.csv(top30,"17cluster-top30.csv")

DoHeatmap(UCCRC,features = top10$gene,assay = "SCT")+NoLegend()
VlnPlot(UCCRC,features = top10$gene[1:10])
VlnPlot(UCCRC,features = c("Creb5","Lcp1","Fpr1"))
table(UCCRC@meta.data$seurat_clusters)
FeaturePlot(UCCRC,features = c("Creb5"),reduction = "tsne")
FeaturePlot(UCCRC,features = c("Fpr1"),reduction = "tsne")
FeaturePlot(UCCRC,features = c("Lcp1"),reduction = "tsne")
UCCRC@meta.data$library=sub("([A-Z]+)\\d+","\\1",UCCRC@meta.data$library)
# 基础版
VlnPlot(
  UCCRC,
  features = "Lcp1",
  group.by = "library",  # 替换为你的分组列名（如"group"）
  pt.size = 0.1,          # 减小点的大小
  cols = c("#1F77B4", "#FF7F0E","red","yellow")  # UC和CRC组的颜色
) + 
  labs(title = "Lcp1 Expression in UC vs CRC") +  # 标题
  theme_classic()         # 简洁主题

# 比较所有细胞中UC vs CRC的表达差异
markers_UC_vs_CRC <- FindMarkers(
  UCandCRC,
  ident.1 = "UC",       # 替换为你的分组标签
  ident.2 = "CRC",
  group.by = "library",
  features = "Lcp1"
)
# 输出结果
markers_UC_vs_CRC
ifnb@meta.data$library=sub("([A-Z]+)\\d+","\\1",ifnb@meta.data$library)
#常规作图
table(ifnb@meta.data$cell_type)
ifnb@meta.data$cell_type=factor(ifnb@meta.data$cell_type,
                                       levels = c("Macrophage",
                                                  "Granulocyte_Neutrophil",
                                                  "Endothelial cell",
                                                  "Enterocyte",
                                                  "CD8 positive alpha-beta T cell",
                                                  "Monocyte",
                                                  "T cell",
                                                  "Goblet cell",
                                                  "B cell",
                                                  "Fibroblast",
                                                  "Pericyte",
                                                  "Plasma cell",
                                                  "Blood vessel endothelial cell",
                                                  "Endothelial cell of lymphatic vessel"))
colour=c("#F08080","#1E90FF","#7CFC00","#808000","#FF00FF","#FA8072","#7B68EE","#800080","#A0522D",
         "#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B")

Idents(ifnb)="cell_type"
markers=FindAllMarkers(object = ifnb,test.use = "wilcox",
                       only.pos = TRUE,
                       logfc.threshold = 0.5,
                       min.pct = 0.5)
all.markers=markers%>%dplyr::select(gene,everything())%>%subset(p_val<0.05)
write.csv(all.markers,"all.markers.csv",row.names = T)
top5=all.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
write.csv(top5,"top5.markers.csv",row.names = T)
top1=all.markers%>%group_by(cluster)%>%top_n(n=1,wt=avg_log2FC)
top5=as.data.frame(top5)
markerdata2=ScaleData(ifnb,features = c("Creb5","Fpr1","Lcp1"),assay = "RNA")

DoHeatmap(markerdata2,
          features = c("Creb5","Fpr1","Lcp1"),
          group.by = "cell_type",
          assay = "RNA",label = F)
VlnPlot(ifnb,features = "Creb5",group.by = "cell_type")
VlnPlot(ifnb,features = "Fpr1",group.by = "cell_type")
VlnPlot(ifnb,features = "Lcp1",group.by = "cell_type")

VlnPlot(ifnb,features = "Creb5",group.by = "cell_type",split.by = "library")
VlnPlot(ifnb,features = "Fpr1",group.by = "cell_type",split.by = "library")
VlnPlot(ifnb,features = "Lcp1",group.by = "cell_type",split.by = "library")

a=c("Creb5","Fpr1","Lcp1")

VlnPlot(
  ifnb,
  features = a,
  pt.size = 0.5,
  group.by = "cell_type",
  stack = TRUE,
  split.by = "library",
  same.y.lims = TRUE,
  flip = T,
  )

DotPlot(ifnb,features = a,group.by = "cell_type")
DotPlot(ifnb,features = a,group.by = "cell_type")+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x = element_text(hjust = 1,vjust = 0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order = 3))+
  theme(axis.text.x = element_text(angle = 90))+
  scale_color_gradientn(colors = c("#330066", "#336699", "#66cc66", "#FFCC33"))


#细胞比例图
table(ifnb$library)
prop.table(table(Idents(ifnb)))
table(Idents(ifnb),ifnb$library)
cellratio=prop.table(table(Idents(ifnb),ifnb$library),margin = 2)
cellratio=as.data.frame(cellratio)
colnames(cellratio)=c("celltype","sample","ratio")

colourCount=length(unique(cellratio$celltype))
ggplot(cellratio)+
  geom_bar(aes(x=sample,y=ratio,fill=celltype),stat = "identity",width = 0.7,size=0.5,colour="#222222")+
  theme_classic()+
  labs(x="Sample",y="Ratio")+
  coord_flip()+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5,linetype = "solid"))+
  scale_fill_manual(values = colour)
#翻转
ggplot(cellratio)+
  geom_bar(aes(x=sample,y=ratio,fill=celltype),stat = "identity",width = 0.7,size=0.5,colour="#222222")+
  theme_classic()+
  labs(x="Sample",y="Ratio")+
 
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5,linetype = "solid"))+
  scale_fill_manual(values = colour)

saveRDS(ifnb,"UCCRC不分组注释后.rds")
