
rm(list=ls(all=TRUE))

library(CIBERSORT)
library(ggplot2)#绘图
library(pheatmap)#绘制热图
library(ggpubr)#堆积比例图
library(reshape2)#数据处理
library(tidyverse)#数据处理


DEG_expr <- read.csv("GSE3629DEG.csv",row.names = 1)
group <- read.csv("3629group.csv",row.names = 1)


boxplot(DEG_expr,outline=F, notch=F , las=2)
qx <- as.numeric(quantile(DEG_expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))  #数据的分布，样本分位数
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
LogC
if (LogC) { 
  DEG_expr[which(DEG_expr <= 0)] <- NaN
  DEG_expr <- log2(DEG_expr) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
}   


library(limma)
DEG_expr=normalizeBetweenArrays(DEG_expr)
boxplot(DEG_expr,outline=FALSE, notch=F , las=2)


LM22_local <- read.table("LM22.txt",header = T,row.names = 1,sep = "\t")
data(LM22)
all(LM22==LM22_local)


result <- cibersort(sig_matrix = LM22, mixture_file = DEG_expr, perm = 1000, QN = TRUE)


result <- as.data.frame(result)
colnames(result)
write.csv(result,"cibersort_result.csv")

result1 <- result[,1:ncol(LM22)]
result1 <- result1[,apply(result1, 2, function(x){sum(x)>0})]

pdf(file="Heatmap.pdf",width = 10,height = 8)
pheatmap(result1,
         color = colorRampPalette(c("#4CB43C", "#FEFCFB", "#ED5467"))(100),
         border="skyblue",
         main = "Heatmap",
         show_rownames = F,
         show_colnames = T,
         cexCol = 1,
         scale = 'row',
         cluster_col=T,
         cluster_row=F,
         angle_col = "90",
         legend = F,
         legend_breaks=c(-3,0,3),
         fontsize_row = 10,
         fontsize_col = 10)
dev.off()



identical(rownames(result1),group$Sample)
data <- cbind(rownames(result1),result1)
colnames(data)[1] <- "Samples"
data <- melt(data,id.vars = c("Samples"))
colnames(data) <- c('Samples','celltype','proportion')

mycolors <- c('#D4E2A7','#88D7A4','#A136A1','#BAE8BC','#C757AF',
              '#DF9FCE','#D5E1F1','#305691','#B6C2E7','#E8EFF7',
              '#9FDFDF','#EEE0F5','#267336','#98CEDD','#CDE2EE',
              '#DAD490','#372E8A','#4C862D','#81D5B0','#BAE8C9',
              '#A7DCE2','#AFDE9C')
pdf(file="stacked bar chart.pdf",width = 10,height = 8)
ggplot(data,
       aes(Samples,proportion,fill=celltype))+geom_bar(stat="identity",position="fill")+#x 轴是变量 Samples，y 轴是变量 proportion，条形的填充颜色由变量 celltype 决定
       scale_fill_manual(values=mycolors)+
       ggtitle("Proportion of immune cells")+theme_gray()+theme(axis.ticks.length=unit(3,'mm'),axis.title.x=element_text(size=11))+
       theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
       guides(fill=guide_legend(title="Types of immune cells"))
dev.off()  



data1 <- cbind(result1,group)
colnames(data1)
data1 <- data1[,c(23,24,1:22)]
rownames(data1) <- NULL
data1 <- pivot_longer(data = data1,
                      cols = 3:22,
                      names_to = "celltype",
                      values_to = "proportion")

pdf(file="group.pdf",width = 10,height = 8)
ggboxplot(data = data1,
          x = "celltype",
          y = "proportion",
          combine = TRUE,
          merge = FALSE,
          color = "black",
          fill = "group",
          palette = NULL,
          title = "TME Cell composition",
          xlab = NULL,
          ylab = "Cell composition",
          bxp.errorbar = FALSE,
          bxp.errorbar.width = 0.2,
          facet.by = NULL,
          panel.labs = NULL,
          short.panel.labs = TRUE,
          linetype = "solid",
          size = NULL,
          width = 0.8,
          notch = FALSE,
          outlier.shape = 20,
          select = NULL,
          remove = NULL,
          order = NULL,
          error.plot = "pointrange",
          label = NULL,
          font.label = list(size = 12, color = "black"),
          label.select = NULL,
          repel = FALSE,
          label.rectangle = TRUE, ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1))  #这个函数的作用是将 x 轴文本旋转 90 度，并调整其对齐方式。
dev.off()  



data2 <- cbind(result1,group)
data2 <- data2[data2$group=='Tumor',]
data2 <- pivot_longer(data = data2,
                      cols = 1:22,
                      names_to = "celltype",
                      values_to = "proportion")
sum(data2$proportion)
pdf(file="one group.pdf",width = 10,height = 8)
ggboxplot(data = data2,
          x = "celltype",
          y = "proportion",
          color = "black",
          xlab = "Types of immune cells",
          ylab = NULL,
          title = "TME Cell composition",
          fill = "celltype",
          legend.position = "bottom",
          ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1)) 
dev.off() 



data3 <- cbind(result1,group)
colnames(data3)
data3 <- data3[,c(23,24,1:22)]
data3 <- pivot_longer(data = data3,
                      cols = 3:24,
                      names_to = "celltype",
                      values_to = "proportion")

pdf(file="group+Pvalue.pdf",width = 10,height = 8)
ggboxplot(data = data3,
          x = "celltype",
          y = "proportion",
          combine = TRUE,
          merge = FALSE,
          color = "black",
          fill = "group",
          palette = c("#ED5432","#81D5B0"),
          title = "TME Cell composition",
          xlab = NULL,
          ylab = "Cell composition",
          bxp.errorbar = FALSE,
          bxp.errorbar.width = 0.2,
          facet.by = NULL,
          panel.labs = NULL,
          short.panel.labs = TRUE,
          linetype = "solid",
          size = NULL,
          width = 0.8,
          notch = FALSE,
          outlier.shape = 20,
          select = NULL,
          remove = NULL,
          order = NULL,
          error.plot = "pointrange",
          label = NULL,
          font.label = list(size = 12, color = "black"),
          label.select = NULL,
          repel = TRUE,
          label.rectangle = TRUE, ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1)) +  #这个函数的作用是将 x 轴文本旋转 90 度，并调整其对齐方式。
          stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.",hide.ns = T) 
dev.off()  
