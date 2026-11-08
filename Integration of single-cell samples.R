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


print(names(alldata@reductions))


alldata.list <- list(UC1,UC2,CRC1,CRC2)


for (i in 1:length(alldata.list)) {
  alldata.list[[i]] <- NormalizeData(alldata.list[[i]], verbose = FALSE)
  alldata.list[[i]] <- FindVariableFeatures(alldata.list[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = FALSE)
}


hvgs_per_dataset <- lapply(alldata.list, function(x) {
  x@assays$SCT@var.features
})
head(hvgs_per_dataset)


venn(hvgs_per_dataset,
     zcolor = scales::hue_pal()(6),  
     opacity = 0.3, ilcs = 1, sncs = 1.2,ilabels = TRUE)

venn::venn(hvgs_per_dataset, opacity = 0.4, zcolor = (scales::hue_pal())(3), cexsn = 1, 
           cexil = 1, lwd = 1, col = "white", frame = F, borders = NA,ilabels = TRUE)


alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:10, 
                                          reduction = "cca")


alldata.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:10, new.assay.name = "CCA")

names(alldata.int@assays)
alldata.int<-SCTransform(alldata.int,vst.flavor = "v2",verbose = T,method="glmGamPoi",vars.to.regress = c("nCount_RNA","percent.mt"))

alldata.int <- RunPCA(alldata.int, npcs = 30, verbose = FALSE)


ElbowPlot_alldata <- ElbowPlot(alldata.int)
ElbowPlot_alldata


par(mar=c(5,5,5,5))
pdf("9.after_integration_PCA_heatmap.pdf",width = 20,height = 20)
DimHeatmap(alldata.int,dims = 1:15,reduction = "pca",cells = 20,ncol = 3)
dev.off()


alldata.int1 <- RunUMAP(alldata.int, dims = 1:8)

DimPlot(alldata.int1,reduction = "umap",pt.size = 1.2)
alldata.int1@meta.data$library=sub("([A-Z]+)\\d+","\\1",alldata.int@meta.data$library)

library(cowplot)

integrating_umap_pca1 <- plot_grid(
  DimPlot(alldata.int1, reduction = "pca", group.by = "library")+NoAxes()+ggtitle("PCA integrated"),
  DimPlot(alldata.int1, reduction = "umap", group.by = "library")+NoAxes()+ggtitle("UMAP integrated")
  
)
integrating_umap_pca1
ggsave("11.pca_umap_integration.pdf",plot = integrating_umap_pca1,device = pdf,width = 25,height = 10,units = "cm",dpi = 300)
library(harmony)

alldata_Harmony <- RunHarmony(alldata.int, "library",max.iter = 200,epsilon.harmony = 1e-4)

names(alldata_Harmony@reductions)

alldata_Harmony <- RunUMAP(alldata_Harmony,  dims = 1:2, 
                           reduction = "harmony")


p1 <- DimPlot(alldata_Harmony,reduction = "umap",label=T )+ggtitle("Harmony")
p2 <- DimPlot(alldata.int1,reduction = "umap",pt.size = 1.2)+ggtitle("UMAP")
p12 <- p1+p2
p12
p1




saveRDS(alldata.int1,"UC+CRC data_after_batcheffect(SCT).rds")
saveRDS(alldata_Harmony,"UC+CRCdata_after_batcheffect(Harmony).rds")
