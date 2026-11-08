


library(limma)
library(ggplot2) 
library(pheatmap) 


expr_data<-read.table("GSE3629.txt",header = T,
                           row.names = 1,sep = "\t")
group<-read.csv("3629.csv",header = T,row.names = 1,sep = ",")


design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr_data)


contrast.matrix <- makeContrasts(Tumor-UC,levels = design)


fit2 <- contrasts.fit(fit, contrast.matrix)  
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1,n = Inf)
DEG$regulate <- ifelse(DEG$adj.P.Val > 0.05, "unchanged",
                   ifelse(DEG$logFC > 1.5, "up-regulated",
                          ifelse(DEG$logFC < -1.5, "down-regulated", "unchanged")))
table(DEG$regulate)
write.table(data.frame(gene_symbol=rownames(DEG),DEG),file = "GSE3629DEG_result.txt",
            sep = "\t",quote = F,row.names = F,col.names = T)


DE_1.5_0.05 <- DEG[DEG$adj.P.Val<0.05&abs(DEG$logFC)>1.5,]
upGene_1.5_0.05 <- DE_1.5_0.05[DE_1.5_0.05$regulate == "up-regulated",]
downGene_1.5_0.05 <- DE_1.5_0.05[DE_1.5_0.05$regulate == "down-regulated",]
write.csv(upGene_1.5_0.05,"upGene(adjp).csv")
write.csv(downGene_1.5_0.05,"downgene(adjp).csv")


pdf("volcano.pdf")
ggplot(DEG,aes(x=logFC,y=-log10(adj.P.Val)))+ 
  geom_point(alpha=0.6,size=3.5,aes(color=regulate))+ 
  ylab("-log10(adj.P.Val)")+ 
  scale_color_manual(values = c("blue", "grey", "red"))+ 
  geom_vline(xintercept = c(-1.5,1.5),lty=4,col ="black",lwd=0.8)+ 
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ 
  theme_bw()  
dev.off()


Genes <- rownames(DEG)
DEG$Genes <- Genes
ggvolcano(data = DEG,x = "logFC",y = "adj.P.Val",label = "Genes",          
          label_number = 10,output = FALSE,
          fills = c("#00AFBB", "#999999", "#FC4E07"),
          colors = c("#00AFBB", "#999999", "#FC4E07"),
          x_lab = "log2FC",
          y_lab = "-Log10 adj.P.Val")



DEG_genes <- DEG[DEG$adj.P.Val<0.05&abs(DEG$logFC)>1.5,]
DEG_gene_expr <- expr_data[rownames(DEG_genes),]
pdf("pheatmap.pdf")
pheatmap(DEG_gene_expr,
         color = colorRampPalette(c("blue","white","red"))(100), 
         scale = "row", 
         border_color = NA, 
         fontsize = 10, 
         show_rownames = F,
         show_colnames = FALSE) 
dev.off()

