k <- Matrix::readMM("matrix.mtx")
library(Seurat)
mtx<-k
dim(mtx)
cl=read.table('barcodes.tsv',
              header = F)[,1]
rl=read.table('features.tsv',
              header = F)[,2]
rownames(mtx) <- rl
colnames(mtx) <- cl
sce=CreateSeuratObject(counts = mtx)
head(cl)

head(rl)

saveRDS(sce,"UC2.rds")



##第二种
library(tidyr)
library(tidyverse)
library(org.Mm.eg.db)
library(scater)
library(scran)
library(pheatmap)
library(limma)
library(edgeR)
library(dynamicTreeCut)
library(cluster)
library(monocle)
aa <- read.csv("GSE122340.csv",row.names = 1,header = T)
all.counts <- as.matrix(aa)
pbmc <- CreateSeuratObject(counts = all.counts, project = "TCGA", names.delim = "-",
                           min.cells = 3, min.features = 200)
