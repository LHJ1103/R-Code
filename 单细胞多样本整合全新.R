library(hdf5r)
library(Seurat) 
library(dplyr)
library(multtest)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(future)
library(glmGamPoi)
library(venn)
UC1[["percent.mt"]]<-PercentageFeatureSet(UC1,pattern = "mt-")
UC2[["percent.mt"]]<-PercentageFeatureSet(UC2,pattern = "mt-")
UC3[["percent.mt"]]<-PercentageFeatureSet(UC3,pattern = "mt-")
CRC1[["percent.mt"]]<-PercentageFeatureSet(CRC1,pattern = "mt-")
CRC2[["percent.mt"]]<-PercentageFeatureSet(CRC2,pattern = "mt-")
CRC3[["percent.mt"]]<-PercentageFeatureSet(CRC3,pattern = "mt-")

#Step3多样本数据整合
#0.加载数据
print(names(alldata@reductions))

#1.分割样本数据构建多个数据集
### 使用SplitObject函数分割seurat对象
alldata.list <- list(UC1,UC2,CRC1,CRC2)

###分别对每个样本数据进行标准化并筛选高变基因
for (i in 1:length(alldata.list)) {
  alldata.list[[i]] <- NormalizeData(alldata.list[[i]], verbose = FALSE)
  alldata.list[[i]] <- FindVariableFeatures(alldata.list[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = FALSE)
}

###提取每个数据集中筛选到的高变基因
hvgs_per_dataset <- lapply(alldata.list, function(x) {
  x@assays$SCT@var.features
})
head(hvgs_per_dataset)

###绘制venn图查看不同样本中高变基因的重叠情况
venn(hvgs_per_dataset,
     zcolor = scales::hue_pal()(6),  # 6种颜色
     opacity = 0.3, ilcs = 1, sncs = 1.2,ilabels = TRUE)

venn::venn(hvgs_per_dataset, opacity = 0.4, zcolor = (scales::hue_pal())(3), cexsn = 1, 
           cexil = 1, lwd = 1, col = "white", frame = F, borders = NA,ilabels = TRUE)

#2.鉴定数据整合的anchors

####使用FindIntegrationAnchors函数鉴定用于多样本数据整合的anchors
alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:10, 
                                          reduction = "cca")

#2.1整合多样本数据集使用IntegrateData函数利用鉴定的anchors进行多样本数据集的整合
alldata.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:10, new.assay.name = "CCA")

names(alldata.int@assays)#查看CCA是否成功添加
alldata.int<-SCTransform(alldata.int,vst.flavor = "v2",verbose = T,method="glmGamPoi",vars.to.regress = c("nCount_RNA","percent.mt"))
#alldata.int <- ScaleData(alldata.int, verbose = FALSE)
alldata.int <- RunPCA(alldata.int, npcs = 30, verbose = FALSE)

##确定数据维度
ElbowPlot_alldata <- ElbowPlot(alldata.int)##胳膊肘确定分群
ElbowPlot_alldata

###PCA热图可视化降维DIM
par(mar=c(5,5,5,5))
pdf("9.after_integration_PCA_heatmap.pdf",width = 20,height = 20)
DimHeatmap(alldata.int,dims = 1:15,reduction = "pca",cells = 20,ncol = 3)
dev.off()##PC=9

#JackStraw图进一步确定分群
#alldata.int <-JackStraw(alldata.int, num.replicate = 100)
#alldata.int <-ScoreJackStraw(alldata.int, dims = 1:4)#结合胳膊肘确定维度dims
#jiackStrawplot1 <- JackStrawPlot(alldata.int, dims = 1:4)#通过保留尽可能少的维度来留存尽可能多的信息

#ggsave(filename = "10.Integration_jiackStrawplot.pdf", plot = jiackStrawplot1, device = 'pdf', width = 27, height = 25, units = 'cm')

####
alldata.int1 <- RunUMAP(alldata.int, dims = 1:8)

DimPlot(alldata.int1,reduction = "umap",pt.size = 1.2)
alldata.int1@meta.data$library=sub("([A-Z]+)\\d+","\\1",alldata.int@meta.data$library)

library(cowplot)
# 可视化
integrating_umap_pca1 <- plot_grid(
  DimPlot(alldata.int1, reduction = "pca", group.by = "library")+NoAxes()+ggtitle("PCA integrated"),
  DimPlot(alldata.int1, reduction = "umap", group.by = "library")+NoAxes()+ggtitle("UMAP integrated")
  #DimPlot(alldata.int1, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE integrated")
)
integrating_umap_pca1
ggsave("11.pca_umap_integration.pdf",plot = integrating_umap_pca1,device = pdf,width = 25,height = 10,units = "cm",dpi = 300)
library(harmony)
####harmony去批次
###Harmony整合
alldata_Harmony <- RunHarmony(alldata.int, "library",max.iter = 200,epsilon.harmony = 1e-4)
#查看Harmony降维是否添加入对象中
names(alldata_Harmony@reductions)
#运行基于样本个数的HAarmony
alldata_Harmony <- RunUMAP(alldata_Harmony,  dims = 1:2, 
                           reduction = "harmony")#依照Harmony运行结果确定降维个数

#可视化
p1 <- DimPlot(alldata_Harmony,reduction = "umap",label=T )+ggtitle("Harmony")
p2 <- DimPlot(alldata.int1,reduction = "umap",pt.size = 1.2)+ggtitle("UMAP")
p12 <- p1+p2
p12
p1
#ggsave("12.Harmonty_remove_batcheffect.pdf",plot = p12,device = pdf,width = 25,height = 10,units = "cm",dpi = 300)



saveRDS(alldata.int1,"UC+CRC（合并）data_after_batcheffect(SCT).rds")
saveRDS(alldata_Harmony,"UC+CRCdata_after_batcheffect(Harmony).rds")
