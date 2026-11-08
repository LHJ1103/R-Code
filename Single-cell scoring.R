library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
inflammatory_gene<-read.csv("geneset.csv")
gene<-as.list(inflammatory_gene)
pbmc<-AddModuleScore(
  object=pbmc,
  features = gene,
  ctrl = 100,
  name="MPs"
)
DimPlot(pbmc)
FeaturePlot(pbmc,features = 'MPs1',reduction = "tsne")

library(ggplot2)
wei<-pbmc@reductions[["umap"]]@cell.embeddings
mydata<-FetchData(pbmc,vars=c("UMAP_1","UMAP_2","MPs1"))
pbmc[["MPs.group"]]<-ifelse(pbmc@meta.data[,"MPs1"]>mean(pbmc@meta.data[,"MPs1"]),"High","Low")
DimPlot(pbmc,group.by = "MPs.group",reduction = "umap")


library(scran)
library(clusterProfiler)
library(org.Mm.eg.db)
View(cc.genes)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

g1_ensembl <- mm.pairs$G1[,1]
g1_symbol <- bitr(g1_ensembl, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
View(g1_symbol)

s_ensembl <- mm.pairs$S[,1]
s_symbol <- bitr(s_ensembl, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
View(s_symbol)

g2m_ensembl <- mm.pairs$G2M[,1]
g2m_symbol <- bitr(g2m_ensembl, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
View(g2m_symbol)


g2m_genes<-g2m_symbol$SYMBOL
g2m_genes<-CaseMatch(search = g2m_genes,match = rownames(pbmc))

s_genes<-s_symbol$SYMBOL
s_genes<-CaseMatch(search = s_genes,match = rownames(pbmc))

g1_genes<-g1_symbol$SYMBOL
g1_genes<-CaseMatch(search = g1_genes,match = rownames(pbmc))


pbmc<-CellCycleScoring(pbmc,g2m.features = g2m_genes,s.features = s_genes,g1.features = g1_genes)
table(pbmc$Phase)
DimPlot(pbmc,group.by = 'Phase',reduction = "tsne")




g1_ensembl <- mm.pairs$G1[,1]
g1_symbol <- bitr(g1_ensembl, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
View(g1_symbol)

s_ensembl <- mm.pairs$S[,1]
s_symbol <- bitr(s_ensembl, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
View(s_symbol)

g2m_ensembl <- mm.pairs$G2M[,1]
g2m_symbol <- bitr(g2m_ensembl, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
View(g2m_symbol)
