setwd("D:/bilibiliR/27、二分类结局变量随机森林筛选变量")

# 下载需要的R包 -----------------------------------------------------------------
# 检查 randomForest 包是否已经安装
if (!require("randomForest", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("randomForest")
}

# 加载所需要的R包 ----------------------------------------------------------------
library(randomForest)#随机森林筛选变量

# 载入示例数据 ------------------------------------------------------------------
mydata <- read.csv("data.csv",row.names = 1)

#简单的数据清洗
mydata <- subset(mydata, select = -time)#删除time这一列
#将二分类变量中的字符转变成数字
colnames(mydata[,1:17])#查看前10列的列名
#因变量
mydata$status <- ifelse(mydata$status == "Tumor",1,0) #死亡为1，存活为0
#自变量转化
mydata$gender <- factor(ifelse(mydata$gender == "female", 0, 1))
mydata$stage <- factor(ifelse(mydata$stage == "early", 0, 1))
mydata$type <- factor(ifelse(mydata$type == "squamous carcinom", 0, 1))
mydata$smoker <- factor(ifelse(mydata$smoker == "nonsmoker", 0, 1))
mydata$alcohol_history <- factor(ifelse(mydata$alcohol_history == "No", 0, 1))

# 定义结局变量和自变量
y <- mydata[, 1]#提取第一列作为二分类结局变量
x <- mydata[, -1]#其余变量为要筛选的自变量


# 随机森林筛选变量 ----------------------------------------------------------------
rf_model <- randomForest(x, # 是一个包含预测变量的矩阵或数据框，用来训练随机森林模型。
                         y, # 是一个向量，包含因变量的观测值，用来训练随机森林模型。
                         ntree = 500, #是一个整数，指定生成的决策树的数量，通常设置为一个较大的值（常用1000），以获得更好的预测性能。
                         importance = TRUE)#是否计算变量的重要性。
# 可视化变量重要性排序结果
varImpPlot(rf_model, 
           type=1, #选择1或者2
           main="Variable Importance",#标题名称
           n.var = 16,#绘制的top变量个数
           scale = T,#是否显示横坐标
           cex = 1.2#字体大小
           )
# 盛夏的果实概念科普，要通俗易懂真的难，看懂点个赞
# MSE（均方差误差），MSE的值越小，说明模型的预测能力越好；
# %IncMSE（Increase in MSE(%)值）:指移除某个特征后，与完整模型相比，平均方差误差（mse）的增加量。
# IncMSE越高，该变量越重要，反之越不重要

varImpPlot(rf_model, 
           type=2, #选择1或者2
           main="Variable Importance",#标题名称
           n.var = 16,#绘制的top变量个数
           scale = T,#是否显示横坐标
           cex = 1.2#字体大小
)
# 盛夏的果实概念科普，要通俗易懂真的难，看懂点个赞
# IncNodePurity:基于基尼系数计算的值,而基尼系数越大,代表分出的类不确定性较大,分类效果不好。
# IncNodePurity用于衡量将变量添加到节点后带来的收益大小，即对节点纯度的提升。
# IncNodePurity越高，该变量越重要，反之越不重要

varImpPlot(rf_model,n.var = 16,cex = 1.2)

# 提取变量重要性
importance <- importance(rf_model)
# ①提取%IncMSE前10的变量
importance <- as.data.frame(importance)
top10_IncMSE <- head(rownames(importance[order(-importance$`%IncMSE`),]), 10)
top10_IncMSE#可以用这个

# ②提取%IncNodePurity前10的变量
top10_ncNodePurity <- head(rownames(importance[order(-importance$IncNodePurity),]), 10)
top10_ncNodePurity#这个更常用一点

# 构建筛选获得变量的表达矩阵
top10_IncMSE_data <- mydata[,top10_IncMSE]
top10_ncNodePurity_data <- mydata[,top10_ncNodePurity]

# 存储结果利于下一步分析 -------------------------------------------------------------
write.csv(top10_IncMSE_data,"top10_IncMSE_data.csv")
write.csv(top10_ncNodePurity_data,"top10_ncNodePurity_data.csv")
