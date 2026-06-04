setwd("D:/bilibiliR/10、提取差异基因的表达谱")

# 加载所需要的R包 ----------------------------------------------------------------
library(limma)
library(ggplot2) #用于绘制火山图
library(pheatmap) #用于绘制热图

# 输入表达矩阵和分组文件 -------------------------------------------------------------
expr_data<-read.table("GSE3629(分析用).txt",header = T,
                      row.names = 1,sep = "\t")
group<-read.csv("3629分组信息.csv",header = T,row.names = 1,sep = ",")

# #构建分组矩阵--design ---------------------------------------------------------
design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr_data)

# #构建比较矩阵——contrast -------------------------------------------------------
contrast.matrix <- makeContrasts(Tumor-UC,levels = design)

# #线性拟合模型构建 ---------------------------------------------------------------
fit <- lmFit(expr_data,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)  
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
DEG <- topTable(fit2, coef = 1,n = Inf)
DEG$regulate <- ifelse(DEG$adj.P.Val > 0.05, "unchanged",
                       ifelse(DEG$logFC > 1.5, "up-regulated",
                              ifelse(DEG$logFC < -1.5, "down-regulated", "unchanged")))
table(DEG$regulate)

# 提取差异基因的表达量 --------------------------------------------------------------
DEG_1_0.05 <- DEG[DEG$adj.P.Val<0.05&abs(DEG$logFC)>1.5,]
deg <- rownames(DEG_1_0.05)
DEG_expr <- expr_data[deg,]
write.csv(DEG_expr,"GSE3629DEG表达谱.csv")
#或者提取任何感兴趣基因的表达量文件
a <- c("ELANE","PADI4","FPR2","IL17A","AZU1","CREB5","CTSG","DYSF","FPR1","LYZ","LCP1","IL1B","MYH9","G0S2","MAPK3","MMP9")
DEG_expr <- expr_data[a,]
write.csv(DEG_expr,"16个DEG表达谱.csv")
