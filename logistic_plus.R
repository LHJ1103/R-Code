setwd("D:/bilibiliR/35、logistic回归列线图美化")

# 下载所需要的R包 ----------------------------------------------------------------
# 检查 rms 包是否已经安装
if (!require("rms", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("rms")
}

# 检查 regplot 包是否已经安装
if (!require("regplot", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("regplot")
}

# 检查 DynNom 包是否已经安装
if (!require("DynNom", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("DynNom")
}

# 加载所需要的R包 ----------------------------------------------------------------
library(rms) # 拟合模型
library(regplot) # 添加个案的列线图
library(DynNom) # 动态列线图

# 载入示例数据 ------------------------------------------------------------------
mydata <- read.csv("data.csv",row.names = 1,sep = ",",header = T)
#变量的处理
mydata$status <- factor(mydata$status,level=c(0,1),labels=c("UC","Tumor"))
mydata$gender <- factor(mydata$gender,level=c(0,1),labels=c("female","male"))
mydata$stage <- factor(mydata$stage,level=c(0,1),labels=c("early","advanced"))

# 开始构建logistic回归模型 ------------------------------------------------------------
#数据打包 
dd = datadist(mydata)
option <- options(datadist = "dd")

#拟合模型
colnames(mydata)#查看列名，选择你要构建模型的变量
formula <- as.formula(status ~ CREB5 + FPR1 + LCP1 )
model <- lrm(formula, # 回归模型的公式,指定自变量和因变量
             data = mydata, # 包含所有变量的数据框
             x=TRUE, # logistic回归也称为"广义线性模型",这个参数指定响应变量的二分类
             y=TRUE) # 参数y也是指定因变量是二分类的

# 结果可视化 -------------------------------------------------------------------
Nomogram_1 <- nomogram(model, #使用nomogram()函数基于Logistic回归模型model生成nomogram对象
                       fun = function(x)1/(1+exp(-x)), #设置nomogram的函数fun为Logistic函数
                       lp=F, #不在nomogram中显示linear predictor的值
                       fun.at = c(0.1,0.5,0.9),#设置nomogram中显示的预测概率值为0.1,0.3,0.5,0.7和0.9 
                       funlabel = "Risk of Tumor")#设置nomogram中预测概率值的标签为“Risk”
#调整一些参数
plot(Nomogram_1,#模型的名称
     xfrac = .30,#变量与图形的占比（调整变量与坐标抽距离）
     cex.var = 1.6,#变量字体加粗
     cex.axis = 1.4,#数轴：字体的大小
     tcl = -0.5,#数轴：刻度的长度
     lmgp = 0.3,#数轴：文字与刻度的距离
     label.every = 1,#数轴：刻度下的文字，1=连续显示，2=隔一个显示一个
     col.grid = c("lightblue","#7CEF7A")#图中线段的颜色
     )

# 添加个案的列线图 ----------------------------------------------------------------
Nomogram_2 <- glm(formula, family=binomial(link = logit), data=mydata)

# 使用 regplot 函数创建一个逻辑回归图，其中包括密度图和箱线图
regplot(Nomogram_2, # Nomogram_2参数指定使用刚才创建的模型
        plots = c("density","boxes") ,# plots 参数指定要绘制的图形类型
        observation = mydata[5,],# 指定在图形中展示的样本，选择第几行
        center = T, # center 参数指定是否将图形的中心点设置为零，对齐变量
        points = TRUE,# points 参数指定是否绘制每个观测值的点
        droplines = F,# droplines 参数指定是否绘制垂直于 x 轴的线
        title = "individual nomogram",# title 参数指定图形的标题
        odds = F, # odds 参数指定是否显示比值或者概率
        interval="confidence",# interval 参数指定使用置信区间还是预测区间
        rank = NULL,# rank 参数指定是否对观测值进行排名
        clickable = F,# 是否启用交互式模式，这个大家一定要换成T，玩一下
        dencol = "skyblue",# dencol 参数指定密度图的颜色
        boxcol = "skyblue"# boxcol 参数指定箱线图的颜色
        )


# 动态列线图 -------------------------------------------------------------------
Nomogram_3 <- glm(formula,data = mydata, family = binomial)

DynNom(Nomogram_3,# Nomogram_3参数指定使用刚才创建的逻辑回归模型
       clevel = 0.55,# clevel 参数指定要使用的置信水平
       DNtitle = "Nomogram",# DNtitle 参数指定图形的标题
       DNxlab = "probability",# DNxlab 参数指定 x 轴的标签
       DNylab = NULL, # DNylab 参数指定 y 轴的标签
       DNlimits = NULL# DNlimits 参数指定 x 轴的限制范围
       )
