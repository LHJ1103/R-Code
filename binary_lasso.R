setwd("D:/bilibiliR/28、二分类结局变量lasso变量筛选")

# 下载需要的R包 -----------------------------------------------------------------
# 检查 glmnet 包是否已经安装
if (!require("glmnet", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("glmnet")
}

# 加载所需要的R包 ----------------------------------------------------------------
library(glmnet)#正则化技术，岭回归、lasso回归和弹性网络

# 载入示例数据 ------------------------------------------------------------------
mydata <- read.csv("data.csv",row.names = 1)

#简单的数据清洗
mydata <- subset(mydata, select = -time)#删除time这一列
#将二分类变量中的字符转变成数字
colnames(mydata[,1:17])#查看前10列的列名
#因变量
mydata$status <- ifelse(mydata$status == "Tumor",1,0) #死亡为1，存活为0
#自变量转化(此处不转化为因子变量)
mydata$gender <- ifelse(mydata$gender == "female", 0, 1)
mydata$stage <- ifelse(mydata$stage == "early", 0, 1)
mydata$type <- ifelse(mydata$type == "squamous carcinom", 0, 1)
mydata$smoker <- ifelse(mydata$smoker == "nonsmoker", 0, 1)
mydata$alcohol_history <- ifelse(mydata$alcohol_history == "No", 0, 1)


# 定义结局变量和自变量(注意区别，此处转化为矩阵)
y <- as.matrix(mydata[, 1])  # 提取第一列作为二分类结局变量
x <- as.matrix(mydata[, 2:17])  # 第2至第51列变量为要筛选的自变量


# Lasso回归筛选变量 -------------------------------------------------------------
set.seed(12345)
lasso_model <- glmnet(x, # 包含自变量的矩阵
                      y, # 因变量向量
                      family = "binomial",# 表示因变量为二元分类变量
                      alpha = 1) # 表示采用L1正则化，即Lasso回归。
print(lasso_model) 
# 结果详解
# Df: 模型的自由度，也就是选择的自变量个数。
# %Dev: 模型的拟合优度，用对数似然比来表示，越接近1说明模型拟合效果越好。
# Lambda: 模型的正则化参数，它是通过交叉验证来选取的，用来控制模型的复杂度，Lambda越大表示模型的正则化强度越大，选择的自变量也越少。

#可视化
plot(lasso_model,# 拟合的Lasso模型
     xvar = "lambda",# 表示将正则化参数lambda作为横坐标
     label = F)# 不在图中显示变量名标签

coef_lasso <- coef(lasso_model,# 拟合的Lasso模型
                   s = 0.071930) # s为Lambda大小，Lambda越大表示模型的正则化强度越大，选择的自变量也越少。
coef_lasso


# 交叉验证选择合适的Lambda ---------------------------------------------------------
cv_model <- cv.glmnet(x, y, family = "binomial",alpha = 1,nfolds = 10)

#可视化
plot(cv_model)

# 根据交叉验证结果，选择lambda值
#我们可以选择平均误差最小的那个λ，即lambda.min
#也可以选择平均误差在一个标准差以内的最大的λ，即lambda.1se。
lambda_min <- cv_model$lambda.min
lambda_min
lambda_1se <- cv_model$lambda.1se
lambda_1se

# 根据lambda值，确定哪些变量应该被保留
coef_cv <- coef(lasso_model, s = lambda_min)
coef_cv 

# 根据回归系数计算OR值
exp(coef_cv)#自然对数的指数函数，它的定义为 exp(x) = e^x，其中 e 是自然对数的基数（约等于 2.718）。


# 把既定Lambda下有意义的变量提取出来 ----------------------------------------------------
coef_cv <- as.matrix(coef_cv)# 先转为矩阵
coef_cv <- data.frame(coef_cv)# 再转为数据框

coef_cv$OR <- exp(coef_cv$s1)# 计算每一个变量的OR值
nonzero_vars <- rownames(coef_cv[coef_cv$OR != 1, ])# OR值不为1的就是我们要筛选的变量
nonzero_vars <- nonzero_vars [2:6]#去除截距项，根据你自己获得的变量数量修改

# 构建筛选获得变量的表达矩阵
lasso_data <- mydata[,nonzero_vars]

# 存储结果利于下一步分析 -------------------------------------------------------------
write.csv(lasso_data,"lasso_data.csv")


#补充内容
# 分割线 ---------------------------------------------------------------------
# 将数据集分成训练集和测试集（适用于样本量较大的情况）
set.seed(12345)
train_index <- sample(1:nrow(mydata), nrow(mydata) * 0.7)# 训练集和测试集7 3分
train_data <- mydata[train_index, ]
test_data <- mydata[-train_index, ]
# 此处为了课程完整性，补充划分数据集的代码
# 示例数据仅152例样本，因此没有划分
# 分割线 ---------------------------------------------------------------------