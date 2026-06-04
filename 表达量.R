setwd("D:/bilibiliR/14、批量绘图")


# 加载所需要的R包 ----------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(Hmisc)#数据处理、分析和可视化/统计分析

# 读取需要处理的数据 ---------------------------------------------------------------
expr_data<-read.csv("16个DEG表达谱.csv",row.names = 1)
group<-read.csv("3629分组信息.csv",header = T,row.names = 1,sep = ",")

# 初步查看所有样本的表达量 ------------------------------------------------------------
mycolors <- c('#D4E2A7','#88D7A4','#A136A1','#BAE8BC','#C757AF',
              '#DF9FCE','#D5E1F1','#305691','#B6C2E7','#E8EFF7',
              '#9FDFDF','#EEE0F5','#267336','#98CEDD','#CDE2EE',
              '#DAD490','#372E8A','#4C862D','#81D5B0','#BAE8C9',
              '#A7DCE2','#AFDE9C')
boxplot(expr_data[,1:59],outline=F, notch=F , las=2,col = mycolors)

# 数据简单处理 ------------------------------------------------------------------
expr_use <- data.frame(t(expr_data)[,1:16])  #转置+提取前10个基因
expr_use <- cbind(group$group,expr_use) #添加分组数据
colnames(expr_use)[1] <- "Group"#修改列名
expr_use <- na.omit(expr_use)#确保没有缺失值
colnames(expr_use)
# 开始绘图 --------------------------------------------------------------------
#箱线图
ggplot(expr_use, aes(x=Group, y=MAPK3,color=Group)) + geom_boxplot() + geom_point(position=position_jitter(width=0.1)) +
       theme(text = element_text(size=12),#调整字体大小和轴标题的大小
             axis.title = element_text(size=18)) + scale_x_discrete(labels=c("Tumor","UC")) + #调整x轴的标签
       theme_minimal(base_size=12, base_family="sans") + 
       scale_color_manual(values = c("#ED545C","#0075C5")) + geom_signif(comparisons=list(c("UC","Tumor")), test = "t.test",map_signif_level=TRUE) 

#一个星号表示P值小于0.05，
#两个星号表示P值小于0.01，
#三个星号表示P值小于0.001。

#小提琴图
ggplot(expr_use, aes(x=Group, y=ELANE,fill=Group)) + geom_violin() + geom_boxplot(outlier.colour = "black",width=0.1,fill="white") +
       theme(text = element_text(size=12),#调整字体大小和轴标题的大小
             axis.title = element_text(size=18)) + 
       scale_x_discrete(labels=c("Tumor","UC")) + #调整x轴的标签
       theme_minimal(base_size=12, base_family="sans") 
       

# 批量绘图 --------------------------------------------------------------------
#定义一个基因列表
gene_list <- colnames(expr_use)[2:17]

#循环建立箱线图
dir.create("expr_compare_boxplot")#创建一个新的文件夹
for (gene in gene_list) {
  print(gene)
  Figure = ggplot(expr_use, aes(x=Group, y=expr_use[,gene],color=Group)) + 
    ylab(gene) + 
    geom_boxplot() +
    geom_point(position=position_jitter(width=0.1)) +
    theme(text = element_text(size=12),#调整字体大小和轴标题的大小
          axis.title = element_text(size=18)) + 
    scale_x_discrete(labels=c("Tumor","UC")) + #调整x轴的标签
    theme_minimal(base_size=12, base_family="sans") + scale_color_manual(values = c("#ED545C","#0075C5")) + geom_signif(comparisons=list(c("UC","Tumor")), test = "wilcox.test",map_signif_level=TRUE)
  ggsave(filename = paste0(gene, ".pdf"), path = "expr_compare_boxplot",plot = Figure,width = 8,height = 6)
}

#循环建立小提琴图
dir.create("expr_compare_violin")#创建一个新的文件夹
for (gene in gene_list) {
  print(gene)
  Figure = ggplot(expr_use, aes(x=Group, y=expr_use[,gene],fill=Group)) + geom_violin() + geom_boxplot(outlier.colour = "black",width=0.1,fill="white") +
    ylab(gene) +
    theme(text = element_text(size=12),#调整字体大小和轴标题的大小
          axis.title = element_text(size=18)) + 
    scale_x_discrete(labels=c("Tumor","UC")) + #调整x轴的标签
    theme_minimal(base_size=12, base_family="sans")
    ggsave(filename = paste0(gene, ".pdf"), path = "expr_compare_violin",plot = Figure,width = 8,height = 6)
}
