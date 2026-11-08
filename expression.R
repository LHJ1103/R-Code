
library(ggplot2)
library(ggpubr)
library(Hmisc)


expr_data<-read.csv("16DEG.csv",row.names = 1)
group<-read.csv("3629group.csv",header = T,row.names = 1,sep = ",")


mycolors <- c('#D4E2A7','#88D7A4','#A136A1','#BAE8BC','#C757AF',
              '#DF9FCE','#D5E1F1','#305691','#B6C2E7','#E8EFF7',
              '#9FDFDF','#EEE0F5','#267336','#98CEDD','#CDE2EE',
              '#DAD490','#372E8A','#4C862D','#81D5B0','#BAE8C9',
              '#A7DCE2','#AFDE9C')
boxplot(expr_data[,1:59],outline=F, notch=F , las=2,col = mycolors)


expr_use <- data.frame(t(expr_data)[,1:16])  
expr_use <- cbind(group$group,expr_use) 
colnames(expr_use)[1] <- "Group"
expr_use <- na.omit(expr_use)
colnames(expr_use)

ggplot(expr_use, aes(x=Group, y=MAPK3,color=Group)) + geom_boxplot() + geom_point(position=position_jitter(width=0.1)) +
       theme(text = element_text(size=12),
             axis.title = element_text(size=18)) + scale_x_discrete(labels=c("Tumor","UC")) + #调整x轴的标签
       theme_minimal(base_size=12, base_family="sans") + 
       scale_color_manual(values = c("#ED545C","#0075C5")) + geom_signif(comparisons=list(c("UC","Tumor")), test = "t.test",map_signif_level=TRUE) 




ggplot(expr_use, aes(x=Group, y=ELANE,fill=Group)) + geom_violin() + geom_boxplot(outlier.colour = "black",width=0.1,fill="white") +
       theme(text = element_text(size=12),
             axis.title = element_text(size=18)) + 
       scale_x_discrete(labels=c("Tumor","UC")) + 
       theme_minimal(base_size=12, base_family="sans") 
       


gene_list <- colnames(expr_use)[2:17]


dir.create("expr_compare_boxplot")
for (gene in gene_list) {
  print(gene)
  Figure = ggplot(expr_use, aes(x=Group, y=expr_use[,gene],color=Group)) + 
    ylab(gene) + 
    geom_boxplot() +
    geom_point(position=position_jitter(width=0.1)) +
    theme(text = element_text(size=12),
          axis.title = element_text(size=18)) + 
    scale_x_discrete(labels=c("Tumor","UC")) + 
    theme_minimal(base_size=12, base_family="sans") + scale_color_manual(values = c("#ED545C","#0075C5")) + geom_signif(comparisons=list(c("UC","Tumor")), test = "wilcox.test",map_signif_level=TRUE)
  ggsave(filename = paste0(gene, ".pdf"), path = "expr_compare_boxplot",plot = Figure,width = 8,height = 6)
}


dir.create("expr_compare_violin")
for (gene in gene_list) {
  print(gene)
  Figure = ggplot(expr_use, aes(x=Group, y=expr_use[,gene],fill=Group)) + geom_violin() + geom_boxplot(outlier.colour = "black",width=0.1,fill="white") +
    ylab(gene) +
    theme(text = element_text(size=12),
          axis.title = element_text(size=18)) + 
    scale_x_discrete(labels=c("Tumor","UC")) + 
    theme_minimal(base_size=12, base_family="sans")
    ggsave(filename = paste0(gene, ".pdf"), path = "expr_compare_violin",plot = Figure,width = 8,height = 6)
}
