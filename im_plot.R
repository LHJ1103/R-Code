# 工作路径 --------------------------------------------------------------------
setwd("D:/bilibiliR/19、基因与免疫细胞相关性")

# 下载所需要的R包 ----------------------------------------------------------------
# 检查 ggcorrplot 包是否已经安装
if (!require("ggcorrplot", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("ggcorrplot")
}

# 检查 tidyr 包是否已经安装
if (!require("tidyr", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("tidyr")
}

# CIBERSORT包安装代码
#install.packages('devtools')
library(devtools)
#if(!require(CIBERSORT))devtools::install_github("Moonerss/CIBERSORT")
#devtools::install_github("Moonerss/CIBERSORT")

# 加载所需要的R包 ----------------------------------------------------------------
library(CIBERSORT)
library(ggcorrplot)
library(tidyr)

# 导入计算cibersort需要的示例数据 ------------------------------------------------------------------
DEG_expr <- read.csv("DEG_expr.csv",row.names = 1)#行名基因名，列名样本名

data <- data(LM22)#cibersort中自带的参考基因集
#假如你运行下面这句代码报错
result <- cibersort(sig_matrix = LM22, #指定用于参考的LM22特征矩阵
                    mixture_file = DEG_expr, #输入待定量的混合物表达矩阵,用于基于LM22矩阵计算细胞类型的丰度
                    perm = 1, # 设置置换次数,用于评估结果的置信区间,默认为1000
                    QN = TRUE # 参数使用快速消元法(Quick Numeric)进行计算,可以显著提高速度,默认为TRUE
                    )
class(DEG_expr)
#解决方式：把输入数据转成矩阵的格式
DEG_expr <- as.matrix(DEG_expr)

# check以下几个地方
# ①输入的DEG_expr文件格式需要是矩阵；
# ②输入的DEG_expr矩阵行名是基因名，列名是样本名；
# ③输入的数据需要是数值型；
# ④输入基因的数目（行名）不能太少；
# ⑤输入的基因（行名）必须是人的基因。cibersort不能用于其他物种的免疫细胞浸润分析。

#cibersort的stat_compare_means函数
# *：表示 p 值小于 0.05，即置信度为 95%。
# **：表示 p 值小于 0.01，即置信度为 99%。
# ***：表示 p 值小于 0.001，即置信度为 99.9%。
# ****：表示 p 值小于 0.0001，即置信度为 99.99%。
# *****：表示 p 值小于 0.00001，即置信度为 99.999%。

# 图形中两组差异不大，检验的P值却很显著，因为这里输入的是免疫细胞浸润占比
# 本身取值较小（同一个样本间所有免疫细胞浸润比值总和为1），样本量（22种免疫细胞）也较小


# 导入上次视频中相关性分析的结果 ---------------------------------------------------------
cor_results <- read.csv("results.csv",row.names = 1,header = T,sep = ",")

head(cor_results)
table(cor_results$Gene)
table(cor_results$im_cell)
# 创建热图
ggplot(cor_results, aes(x = Gene, y = im_cell, fill = Cor)) +   
  # 使用ggplot初始化绘图,设置映射关系,# x轴为Gene,y轴为im_cell,填充颜色映射到Cor相关系数
  geom_tile() +   
  # 添加热力图平铺图层
  scale_fill_gradient2(low = "#7AA7EF", mid = "white", high = "#E9384E") +
  # 设置二段渐变填充颜色,从蓝到白再到红
  labs(title = "Correlation between immune cells and genes", x = "Gene", y = "Immune Cell") +  
  # 添加标题和坐标轴标签
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  # 设置x轴标签textelement属性,包括旋转90度,垂直居中等
  geom_text(aes(label = ifelse(p.value < 0.05, format(p.value, nsmall = 2), "")), size = 3, color = "white", hjust = 0.5) +   
  # 添加p值文本,当p值<0.05时显示p值,否则显示空字符串,设置文本大小、颜色和水平对齐方式
  scale_y_discrete(limits = rev(levels(cor_results$im_cell))) 
  # 设置y轴为离散型,并倒序显示im_cell级别


# 创建热图2（添加*版） -------------------------------------------------------------
ggplot(cor_results, aes(x = Gene, y = im_cell, fill = Cor)) +
  geom_tile() +
  scale_fill_gradient2(low = "#7AA7EF", mid = "white", high = "#E9384E") +
  labs(title = "Correlation between immune cells and genes", x = "Gene", y = "Immune Cell") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # 根据p值使用不同的符号
  geom_text(aes(label = ifelse(p.value < 0.001, "***",  
                               ifelse(p.value < 0.01, "**",
                                      ifelse(p.value < 0.05, "*","")))), size = 3, color = "white", hjust = 0.5) +
  scale_y_discrete(limits = rev(levels(cor_results$im_cell))) +
  # 对图中的字体进行美化、修改
  theme(axis.text.x = element_text(size = 15, face = "bold"),  
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title=element_text(size=20,face="bold"),
        panel.background = element_rect(fill = "white", color = "grey"),
        strip.background = element_rect(fill = "lightgray", color = "grey"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey", size = 0.5),
        panel.margin = unit(5, "lines"))  

  
  

  