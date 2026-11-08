
library(limma)
library(ggplot2) 
library(pheatmap) 


expr_data<-read.table("GSE3629.txt",header = T,
                      row.names = 1,sep = "\t")
group<-read.csv("3629group.csv",header = T,row.names = 1,sep = ",")


design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr_data)


contrast.matrix <- makeContrasts(Tumor-UC,levels = design)


fit <- lmFit(expr_data,design) 
fit2 <- contrasts.fit(fit, contrast.matrix)  
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1,n = Inf)
DEG$regulate <- ifelse(DEG$adj.P.Val > 0.05, "unchanged",
                       ifelse(DEG$logFC > 1.5, "up-regulated",
                              ifelse(DEG$logFC < -1.5, "down-regulated", "unchanged")))
table(DEG$regulate)


DEG_1_0.05 <- DEG[DEG$adj.P.Val<0.05&abs(DEG$logFC)>1.5,]
deg <- rownames(DEG_1_0.05)
DEG_expr <- expr_data[deg,]
write.csv(DEG_expr,"GSE3629DEG.csv")

a <- c("ELANE","PADI4","FPR2","IL17A","AZU1","CREB5","CTSG","DYSF","FPR1","LYZ","LCP1","IL1B","MYH9","G0S2","MAPK3","MMP9")
DEG_expr <- expr_data[a,]
write.csv(DEG_expr,"16个DEG.csv")
