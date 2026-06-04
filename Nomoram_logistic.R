# 加载工作路径 ------------------------------------------------------------------
setwd("D:/bilibiliR/30、Nomogram_logistic")

# 下载所需要的R包 ----------------------------------------------------------------
# 检查 VIM 包是否已经安装
if (!require("VIM", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("VIM")
}

# 检查 rms 包是否已经安装
if (!require("rms", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("rms")
}

# 检查 nomogramFormula 包是否已经安装
if (!require("nomogramFormula", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("nomogramFormula")
}

# 检查 pROC 包是否已经安装
if (!require("pROC", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("pROC")
}

# 检查 rmda 包是否已经安装
if (!require("rmda", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("rmda")
}

# 加载需要的R包 -----------------------------------------------------------------
library(VIM)#缺失值可视化
library(rms)#拟合模型
library(nomogramFormula)#计算列线图得分
library(pROC)#绘制ROC曲线，计算AUC和95%置信区间
library(rmda)#临床决策曲线和临床影响曲线
library(readr)
library(foreign)
library(rms)
library(tableone)
library(broom)
library(regplot)
library(ResourceSelection)
library(nricens)
library(PredictABEL)
library(ggDCA)
library(caret)
library(survivalROC)

# 载入示例数据 ------------------------------------------------------------------
mydata <- read.csv("data.csv",row.names = 1,sep = ",",header = T)

val.data = read.csv("GSE37283YZ.csv",row.names = 1,sep = ",",header = T)
val.data2 = read.csv("E:/p小于0.05 Fc大于1.5/列线图/GSE206285YZ.csv",row.names = 1,sep = ",",header = T)
val.data3 = read.csv("E:/p小于0.05 Fc大于1.5/列线图/GSE75214YZ.csv",row.names = 1,sep = ",",header = T)
val.data4 = read.csv("E:/p小于0.05 Fc大于1.5/列线图/GSE59071YZ.csv",row.names = 1,sep = ",",header = T)
mydata_copy <- mydata#数据备份
# 查看数据和简单数据清洗 -------------------------------------------------------------
# 查看数据中的缺失值
aggr(mydata, col=c('skyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(mydata), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Percentage"))
# 删除缺失值所在的行
mydata <- na.omit(mydata)
summary(mydata)
#生存变量的处理
table(mydata$status)#检查生存状态
mydata$status <- factor(mydata$status,level=c(0,1),labels=c("UC","Tumor"))
#最终变量不宜过多，一般 变量数＜阳性时间/10-20

#二分类变量的转换（转换成因子，列线图能展示相应分组字符，计算评分时不要运行）
mydata$gender <- factor(mydata$gender,level=c(0,1),labels=c("female","male"))
mydata$stage <- factor(mydata$stage,level=c(0,1),labels=c("early","advanced"))

# 开始构建logistic回归模型 ------------------------------------------------------------
#数据打包 
dd = datadist(mydata)
option <- options(datadist = "dd")
model.all = glm(status ~ CREB5 + FPR1 + LCP1 ,
                data = mydata,
                family = binomial(link = "logit"))
summary(model.all)
exp(cbind("OR"=coef(model.all),confint(model.all)))
model.all = glm(status ~ CREB5 + FPR1 + LCP1 ,
                data = mydata,
                family = binomial(link = "logit"),
                x = TRUE)
ddist = datadist(mydata)
options(datadist = "ddist")
model.all2 = lrm(status ~ CREB5 + FPR1 + LCP1 ,
                data = mydata,
                x = TRUE,
                y = TRUE)
Nomogram_1 <- nomogram(model.all2,
                       fun = function(x)1/(1+exp(-x)),
                       lp=F,
                       fun.at = c(0.1,0.3,0.5,0.7,0.9),
                       funlabel = "Risk of UC")

plot(Nomogram_1)
model.all = glm(status ~ CREB5 + FPR1 + LCP1 ,
                data = mydata,
                family = binomial(link = "logit"))
regplot(model.all,
        plots = c("density","boxes"),
        observation = mydata[12,],
        clickable = FALSE,
        title = "individual nomogram",
        dencol = "skyblue",# dencol 参数指定密度图的颜色
        boxcol = "skyblue")# boxcol 参数指定箱线图的颜色)
mydata$pred.all = predict(object = model.all,
                          newdata = mydata,
                          type = "response")

val.data$pred.all = predict(object = model.all,
                            newdata = val.data,
                            type = "response")


val.data2$pred.all = predict(object = model.all,
                            newdata = val.data2,
                            type = "response")


val.data3$pred.all = predict(object = model.all,
                            newdata = val.data3,
                            type = "response")

val.data4$pred.all = predict(object = model.all,
                             newdata = val.data4,
                             type = "response")

val.roc = roc(status~pred.all,data = val.data,smooth = FALSE)
plot(val.roc,
     print.auc = TRUE,
     print.thres = TRUE,
     main = "ROC CURVE OF GSE37283",
     col = "#493D26",
     print.thres.col = "#493D26",
     identity.lty = 2,
     identity.lwd = 2)

ROC = roc(val.data2$status,val.data2$pred.all)
auc(ROC)
ci(auc(ROC))
plot(1-ROC$specificities,  
     #x轴为1-Specificity,特异度的补值
     ROC$sensitivities,type = "l", 
     #y轴为ROC$sensitivities,敏感度
     col="red", 
     #线条颜色为红色
     lty=1,
     #线型为实线
     xlab = "1-Specificity",
     #x轴标签为1-Specificity,特异度的补值
     ylab = "Sensitivity",
     #y轴标签为Sensitivity,敏感度
     lwd=3) 
#线条宽度为2
abline(0,1) 
#添加对角线
legend(0.3,0.1,
       c("AUC of GSE206285 : 0.9595"),
       #图例位置及内容
       lty = c(1),   
       #图例线型为实线
       lwd = c(3),   
       #图例线条宽度为2
       col = c("red"), 
       #图例颜色为红色
       bty = "0")    


val.roc2 = roc(status~pred.all,data = val.data2,smooth = FALSE)
plot(val.roc2,
     print.auc = TRUE,
     print.thres = TRUE,
     main = "ROC CURVE OF GSE206285",
     col = "red",
     print.thres.col = "red",
     identity.lty = 2,
     identity.lwd = 1)

val.roc3 = roc(status~pred.all,data = val.data3,smooth = FALSE)
plot(val.roc3,
     print.auc = TRUE,
     print.thres = TRUE,
     main = "ROC CURVE OF GSE75214",
     col = "blue",
     print.thres.col = "blue",
     identity.lty = 2,
     identity.lwd = 2)

ROC4 = roc(val.data4$status,val.data4$pred.all)
auc(ROC4)
ci(auc(ROC4))
plot(1-ROC4$specificities,  
     #x轴为1-Specificity,特异度的补值
     ROC4$sensitivities,type = "l", 
     #y轴为ROC$sensitivities,敏感度
     col="red", 
     #线条颜色为红色
     lty=1,
     #线型为实线
     xlab = "1-Specificity",
     #x轴标签为1-Specificity,特异度的补值
     ylab = "Sensitivity",
     #y轴标签为Sensitivity,敏感度
     lwd=3) 
#线条宽度为2
abline(0,1) 
#添加对角线
legend(0.3,0.1,
       c("AUC of GSE59071 : 0.9597"),
       #图例位置及内容
       lty = c(1),   
       #图例线型为实线
       lwd = c(3),   
       #图例线条宽度为2
       col = c("red"), 
       #图例颜色为红色
       bty = "0")    


   

val.roc4 = roc(status~pred.all,data = val.data4,smooth = FALSE)
plot(val.roc4,
     print.auc = TRUE,
     print.thres = TRUE,
     main = "ROC CURVE OF GSE59071",
     col = "blue",
     print.thres.col = "blue",
     identity.lty = 2,
     identity.lwd = 2)

plot(val.roc,
     print.auc = TRUE,
     print.thres = TRUE,
     main = "ROC CURVE",
     col = "blue",
     print.thres.col = "blue",
     identity.lty = 2,
     identity.lwd = 2)
plot(val.roc2, add=TRUE, col="#99DD98",print.auc = TRUE,print.thres = TRUE,print.thres.col = "#99DD98",identity.lty = 2,identity.lwd = 2)
plot(val.roc4, add=TRUE, col="#8F94E0")
plot(roc4, add=TRUE, COL="#566D7E")
plot(roc5, add=TRUE, COL="#50EBEC")
plot(roc6, add=TRUE, COL="#FFFF00")
plot(roc7, add=TRUE, COL="#493D26")
plot(roc8, add=TRUE, COL="#C5908E")
plot(roc9, add=TRUE, COL="#6A287E")
plot(roc10, add=TRUE, COL="#7D0552")

#拟合模型
colnames(mydata)#查看列名，选择你要构建模型的变量
formula <- as.formula(status ~ CREB5 + FPR1 + LCP1 )
model <- lrm(formula, # 回归模型的公式,指定自变量和因变量
           data = mydata, # 包含所有变量的数据框
           x=TRUE, # logistic回归也称为"广义线性模型",这个参数指定响应变量的二分类
           y=TRUE) # 参数y也是指定因变量是二分类的
model#查看模型具体情况
#1. 模型概况:显示数据集obs数(152)、LR卡方值(29.32)和p值(0.0001),说明模型整体是显著的。R2值为0.239,表现一定的预测能力。
#2. 区分度指数:包括R2(决定系数)、Brier分数(预测精度)、gamma(区分度)、tau-a(区分度)等,用于评估模型的预测能力。值越大越好。
#3. 变量系数:每一列代表一个自变量,从左至右是:
#- Coef:变量的Logistic回归系数,代表变量的权重大小。正值意味着当变量增加时,1的对数几率增加,负值相反。
#- S.E.:系数的标准误,用于后续的Wald检验。
#- Wald Z:系数显著性的Wald统计量,Z越大绝对值越显著。
#- Pr(>|Z|):对应Z值的p值。小于0.05为显著变量。

OR <- exp(model$coefficients)#计算Logistic回归模型中每个自变量的比值比(Odds Ratio, OR)
OR

# 结果可视化 -------------------------------------------------------------------
Nomogram_1 <- nomogram(model,
                       fun = function(x)1/(1+exp(-x)),
                       lp=F,
                       fun.at = c(0.1,0.5,0.9),
                       funlabel = "Risk of Tumor")
plot(Nomogram_1)
# 答疑：
# 在Nomogram图的绘制中,线段长度同时取决于两个因素:
# 1) 变量的回归系数:系数绝对值越大的变量,其线段越长,其他条件相同。
# 2) 变量的取值范围:范围越大的变量,其线段也越长,其他条件相同。

#调整一些参数
plot(Nomogram_1,
     #模型的名称
     xfrac = .35,
     #变量与图形的占比（调整变量与坐标抽距离）
     cex.var = 1,
     #变量字体加粗
     cex.axis = 1,
     #数轴：字体的大小
     tcl = -0.5,
     #数轴：刻度的长度
     lmgp = 0.3,
     #数轴：文字与刻度的距离
     label.every = 1,
     #数轴：刻度下的文字，1=连续显示，2=隔一个显示一个
     col.grid = gray(c(0.8,0.95)))





#计算列线图评分(要计算评分，数据中不能含有因子变量)
options(option) 
#设置options()函数相关参数
results <- formula_rd(nomogram = Nomogram_1)
#使用formula_rd()函数基于nomogram对象Nomogram_1生成公式formula和相应的数据框rd
val.data$points<-points_cal(formula = results$formula,rd=val.data)
#使用points_cal()函数基于formula公式和原始数据框mydata计算points列,并添加到mydata
head(mydata) 
#查看mydata的数据框的前6行

#可以根据列线图评分（中位数、均值、四分位数）对所有样本进行分组

# 模型评估 -----------------------------------------------------------------
# ROC曲线 -------------------------------------------------------------------
#AUC计算
install.packages("brglm")  # 如果未安装
library(brglm)

# 替换为你的 formula 和 mydata
model_ROC <- brglm(formula, data = mydata, family = binomial())
summary(model_ROC)

model_ROC <- glm(formula,data = mydata,family = binomial())
#lrm()函数属于Design包,专门用于拟合Logistic回归模型。
#glm()函数属于stats包,是用于拟合广义线性模型(Generalized Linear Models)的泛用函数,可以拟合Logistic回归、Poisson回归等多种模型。

#type=“response”给出具体的预测概率，而type=“class”按规定的阈值给出分类
mydata$predvalue <- predict(model_ROC,type="response")
# 计算AUC和95%CI
ROC <- roc(mydata$status,mydata$predvalue)
auc(ROC)
ci(auc(ROC))
#绘制ROC曲线
plot(1-ROC$specificities,  
     #x轴为1-Specificity,特异度的补值
     ROC$sensitivities,type = "l", 
     #y轴为ROC$sensitivities,敏感度
     col="red", 
     #线条颜色为红色
     lty=1,
     #线型为实线
     xlab = "1-Specificity",
     #x轴标签为1-Specificity,特异度的补值
     ylab = "Sensitivity",
     #y轴标签为Sensitivity,敏感度
     lwd=3) 
     #线条宽度为2
abline(0,1) 
#添加对角线
legend(0.5,0.25,
       c("AUC of Nomogram:0.8808"),
       #图例位置及内容
       lty = c(1),   
       #图例线型为实线
       lwd = c(3),   
       #图例线条宽度为2
       col = c("red"), 
       #图例颜色为红色
       bty = "0")    

# 分割线 ---------------------------------------------------------------------
#假如你想绘制其中几个变量的ROC曲线，并在同一张图中展示
colnames(mydata)
roc1 <- roc(mydata$status,mydata$predvalue)#模型的ROC
roc2 <- roc(mydata$status,mydata$CREB5)#gender变量的ROC
roc3 <- roc(mydata$status,mydata$FPR1)#KCNJ13变量ROC
roc4 <- roc(mydata$status,mydata$LCP1)
roc5 <- roc(mydata$status,mydata$S100A8)
plot(roc1, 
     #roc1为ROC对象,用于绘制ROC曲线图
     col="red",  
     #线条颜色为红色
     print.auc=F, 
     #在图中打印AUC值
     print.thres=F,
     #在图中打印最佳切割点的阈值
     lty=1,  
     #线型为实线
     xlab = "Specificity",
     #x轴标签为Specificity,特异度
     ylab = "Sensitivity", 
     #y轴标签为Sensitivity,敏感度
     lwd=3)  
     #线条宽度为2
plot(roc2, add=TRUE, col="#99DD98",lwd=2.5)
plot(roc3, add=TRUE, col="#8F94E0",lwd=2.5)
plot(roc4, add=TRUE, col="#6A287E",lwd=2.5)
plot(roc5, add=TRUE, COL="#50EBEC")
plot(roc6, add=TRUE, COL="#FFFF00")
plot(roc7, add=TRUE, COL="#493D26")
plot(roc8, add=TRUE, COL="#C5908E")
plot(roc9, add=TRUE, COL="#6A287E")
plot(roc10, add=TRUE, COL="#7D0552")
#简单计算另外两个变量的AUC和95%CI
auc(roc1)
ci(auc(roc1))

auc(roc2)
ci(auc(roc2))
auc(roc3)
ci(auc(roc3))
auc(roc4)
ci(auc(roc4))

#添加图例
legend(0.72,0.28,
       legend=c("Model:0.8808(95% CI: 0.7711-0.9905)",
                "CREB5:0.7558(95% CI: 0.5978-0.9138)",
                "FPR1:0.8169(95% CI: 0.6881-0.9456)",
                "LCP1:0.8343(95% CI: 0.7084-0.9602)"),
       col=c("red","#99DD98","#8F94E0","#6A287E"),
       lty=1,
       lwd=2,
       cex = 1)

# 校准曲线 ---------------------------------------------------------------------
cal <- calibrate(model,method ="boot" ,B=1000)
# 使用calibrate()函数基于bootstrap方法对模型model进行校准,bootstrap重复次数为1000

# 设置图形参数
par(cex.lab = 1.2, cex.axis = 1.2, lwd = 2)

# 绘制校准曲线
plot(cal, 
     xlim = c(0, 1),  
     ylim = c(0, 1), 
     xlab = "Predicted Probability",
     ylab = "Observed Probability", 
     subtitles = FALSE)

# 添加参考线
abline(0, 1, col = "red", lty = 4)

# 添加图例
legend(0.001,1.1, 
       legend = c("Calibration Curve", "Ideal Line"), 
       col = c("black", "red"), 
       lty = c(1, 2), 
       lwd = 2)




plot(cal, 
     #绘制cal校准对象的校准曲线图
     xlim=c(0,1),  
     #x轴范围从0到1
     ylim=c(0,1), 
     #y轴范围从0到1
     xlab="Predicted Probability",
     #x轴标签为Predicted Probability,预测概率
     ylab="Observed Probability", 
     #y轴标签为Observed Probability,观察概率
     subtitles=FALSE) 
     #不显示子标题



# 临床决策曲线和临床影响曲线 -----------------------------------------------------------
model_DCA <- decision_curve(formula, 
                            #公式formula为模型的公式 
                            data = mydata_copy, 
                            #data为模型拟合的数据集 
                            family = binomial(link = 'logit'), 
                            #模型属于二项式Logistic回归  
                            thresholds = seq(0,1,by=0.01),  
                            #阈值范围从0到1,间隔为0.01
                            confidence.intervals = 0.95, 
                            #置信区间为95%
                            study.design = 'case-control', 
                            #研究设计为病例对照研究 
                            population.prevalence = 0.3)  
                            #人群基线发病概率为0.3
#DCA曲线
plot_decision_curve(model_DCA,  
                    #绘制model_DCA的决策曲线图 
                    curve.names = c('Nomogram'), 
                    #曲线名称为Nomogram 
                    xlim = c(0,1.0),   
                    #x轴范围0到1 
                    cost.benefit.axis = FALSE,  
                    #不显示成本效益坐标轴 
                    col = c('red'),   
                    #曲线颜色为红色
                    confidence.intervals = FALSE, 
                    #不显示置信区间 
                    standardize = FALSE) 
                    #不进行标准化
#临床影响曲线
plot_clinical_impact(model_DCA, 
                     #绘制临床影响曲线图
                     population.size = 1000,  
                     #总人群规模为1000人 
                     cost.benefit.axis = T,    
                     #显示成本效益坐标轴 
                     n.cost.benefits = 8, 
                     #成本效益轴分为8段 
                     col = c('red','blue'), 
                     #红色代表病例,蓝色代表对照 
                     confidence.intervals = F) 
                     #不显示置信区间