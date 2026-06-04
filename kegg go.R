setwd("E:/Rdata/KEGG GO/")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

upGene = read.csv("16个DEG表达谱.csv",row.names = 1)
#Genes = read.csv("E:/Rdata/KEGG GO/Genes.csv")
Genes = bitr(rownames(upGene),
            fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
#write.csv(Genes,"Genes.csv")

GO = enrichGO(gene = Genes$ENTREZID,
              OrgDb = org.Hs.eg.db,
              keyType = "ENTREZID",
              ont = "all",
              pAdjustMethod = "BH",
              pvalueCutoff = 1,
              qvalueCutoff = 1,
              minGSSize = 5,
              maxGSSize = 5000,
              readable = TRUE)
write.table(data.frame(ID=row.names(GO@result),GO@result),
            file = "GO_enrich_result.txt",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
pdf(file = "GO柱状图.pdf",width = 10,height = 8)
barplot(GO,drop = TRUE,
        showCategory = 6,split = "ONTOLOGY") + facet_grid(ONTOLOGY~.,scales = "free")
dev.off()  

KEGG = enrichKEGG(gene = Genes$ENTREZID,
                  organism = "hsa",
                  keyType = "kegg",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)  
barplot(KEGG,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 15,
        title = "KEGG_enrichment")
dotplot(KEGG)
View(as.data.frame(KEGG))
#browseKEGG(KEGG,"hsaxxxxxx")

KEGG_results = DOSE::setReadable(KEGG,
                                 OrgDb = org.Hs.eg.db,
                                 keyType = "ENTREZID")
write.table(data.frame(ID=rownames(KEGG_results@result),KEGG_results@result),
            file = "KEGG_enrich_result.txt",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
