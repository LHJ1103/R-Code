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

# 加载需要的R包 -----------------------------------------------------------------
library(ggcorrplot)
library(tidyr)

# 导入数据 --------------------------------------------------------------------
data <- read.csv("data.csv",row.names = 1)#行名基因名，列名样本名
# cibersort运算结果
cibersort_result <- read.csv("cibersort_result.csv",row.names = 1)#行名样本名，列名免疫细胞浸润评分
#数据格式整理
data <- t(DEG_expr)#数据转置，行名是样本名，列名基因名
cibersort_result <- cibersort_result[,1:22]

# 解决历史遗留问题 ----------------------------------------------------------------
#问题一：对于Cibersort免疫细胞浸润分析，如果输入的基因表达矩阵是以log转换的形式表示的，
#需要进行反转换以获得原始表达值，因为Cibersort需要原始的表达值来进行细胞比例估计。

#进行log逆转换(此处以2为底举例，实际以GEO官网数据集介绍为准)
if(max(data) < 50) {data <- 2^data}

#问题二：在Cibersort中，输入基因表达矩阵应该是进行标准化处理后的数据。
#如果对数转换后再进行反转换，得到的数据是未标准化的原始数据，这与Cibersort需要的标准化数据是不一样的。
data <- t(DEG_expr)#数据转置，行名是样本名，列名基因名
#check一下
identical(rownames(data),rownames(cibersort_result))

#数据简单处理
data <- data[,1:3]#提取前10个基因
data <- as.data.frame(data)
# Shapiro-Wilk检验
shapiro.test(data$SHISA3)#适用于小样本数据，N≤50
# Kolmogorov-Smirnov检验(柯尔莫哥洛夫检验)
ks.test(data$LCP1, pnorm, mean = mean(data$LCP1), sd = sd(data$LCP1))#适用于大样本数据，N＞50

# 计算相关性 -------------------------------------------------------------------
colnames(data)
colnames(cibersort_result)
target_Gene <- "CREB5"#确定目标基因
expr <- as.numeric(data[,target_Gene])

cor.test(cibersort_result$T.cells.CD8, expr, use = "everything", method = "spearman")

# 封装一个函数来快速计算 -------------------------------------------------------------
# spearman相关性分析
calculate_correlation <- function(Gene_expr, cibersort_result) {
  cor_matrix <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:ncol(cibersort_result)) {
    result <- cor.test(Gene_expr, cibersort_result[, i], method = "spearman")
    new_row <- data.frame("Gene" = "Gene", "im_cell" = colnames(cibersort_result)[i], "Cor" = result$estimate, "p-value" = result$p.value)
    cor_matrix <- rbind(cor_matrix, new_row)
  }
  return(cor_matrix)
}
# pearson相关性分析
calculate_correlation <- function(Gene_expr, cibersort_result) {
  cor_matrix <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:ncol(cibersort_result)) {
    result <- cor.test(Gene_expr, cibersort_result[, i], method = "pearson")
    new_row <- data.frame("Gene" = "Gene", "im_cell" = colnames(cibersort_result)[i], "Cor" = result$estimate, "p-value" = result$p.value)
    cor_matrix <- rbind(cor_matrix, new_row)
  }
  return(cor_matrix)
}

# 计算5个基因与免疫细胞浸润分数的相关性 -----------------------------------------------------
# 选择要计算相关性的列
data_to_calculate <- data[,1:3]
# 新建一个空的数据框保存结果
results <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
# 使用for循环遍历数据框中的每一列，并计算相关性
for (i in 1:ncol(data_to_calculate)) {
  print(i)
  gene_expr <- data_to_calculate[, i]
  corr_result <- calculate_correlation(gene_expr, cibersort_result)
  # 将每次计算的结果添加到新的数据框中
  results <- rbind(results, corr_result)
}
# 查看列名，下方修改的依据
colnames(data_to_calculate)
# 手动修改每一个基因名称
results$Gene <- c(rep("CREB5", 22),
                  rep("FPR1", 22),
                  rep("LCP1", 22))
                  
write.csv(results,"results.csv")
results <- results[,c(1,2,3)]#转换数据，这里不需要显著性p值，也可以通过显著性P值进行筛选
colnames(results)
# 使用spread函数将长数据转换为宽格式
results_wide <- spread(results, key = im_cell, value = Cor)
# 将基因名设置为行名
rownames(results_wide) <- results_wide$Gene
# 删除原始数据框中的基因名列
results_wide$Gene <- NULL

# 绘制相关性矩阵图
ggcorrplot(t(results_wide), 
           # 显示颜色图例
           show.legend = T, 
           # 设置相关系数矩阵的颜色，其中第一个为最小值对应的颜色，
           # 中间的白色为0对应的颜色，最后一个为最大值对应的颜色
           colors = c("#2166AC", "white", "#B2182B"), 
           # 设置数字显示的位数
           digits = 2, 
           # 显示变量标签（默认为TRUE）
           lab = T)
           
  


# 另一种简单的方法 （感受R包的魅力）----------------------------------------------
# 计算相关系数矩阵
corr_mat <- cor(data, cibersort_result, method = "spearman")
# 绘制相关性矩阵图
ggcorrplot(corr_mat, 
           # 显示颜色图例
           show.legend = T, 
           # 设置相关系数矩阵的颜色，其中第一个为最小值对应的颜色，
           # 中间的白色为0对应的颜色，最后一个为最大值对应的颜色
           colors = c("#2166AC", "white", "#B2182B"), 
           # 设置数字显示的位数
           digits = 2, 
           # 显示变量标签（默认为TRUE）
           lab = T) + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#修改横轴标签方向
