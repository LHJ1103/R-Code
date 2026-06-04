# 安装R包
library(xgboost) #用于拟合xgboost模型
library(caret)   #用于一般的数据准备和模型拟合
library(pROC) #制作ROC曲线
#加载数据
train <- read.csv("data - 副本.csv",row.names = 1)
#查看数据结构
str(train)
#让这个例子可重现
set.seed(1234)
#分为训练集(60%)和测试集(40%)
parts = createDataPartition(train$fen.lasso, p = .6, list = F)
train1 = train[parts, ]
test1 = train[-parts, ]
#在训练集中定义预测变量和响应变量
train_x = data.matrix(train1[, -1])
train_y = train1[,1]
#在测试集中定义预测变量和响应变量
test_x = data.matrix(test1[, -1])
test_y = test1[, 1]
#定义最终的训练集和测试集
xgb_train = xgb.DMatrix(data = train_x, label = train_y)
xgb_test = xgb.DMatrix(data = test_x, label = test_y)

#定义观察名单中
watchlist = list(train=xgb_train, test=xgb_test)
#拟合XGBoost模型并在每轮显示训练和测试数据
model = xgb.train(data = xgb_train, max.depth = 3, watchlist=watchlist, nrounds = 100)


#定义最终模型
final = xgboost(data = xgb_train, max.depth = 3, nrounds = 7, verbose = 0)
library("DALEX")

predict_logit <- function(model, x) {
  raw_x <- predict(model, x)
  exp(raw_x)/(1 + exp(raw_x))
}

logit <- function(x) {
  exp(x)/(1+exp(x))
}

explainer_xgb <- explain(final,
                         data = train_x,
                         y = as.numeric(train_y),
                         predict_function = predict_logit,
                         link = logit,
                         label = "xgboost")
explainer_xgb

sv_xgb_satisfaction_level <-DALEX::model_profile(explainer_xgb,
                                                 type = "partial")
head(sv_xgb_satisfaction_level)
plot(sv_xgb_satisfaction_level)

library("breakDown")
library(ggplot2)
nobs <- train_x[1L, , drop = FALSE]

explain_2 <- broken(final, new_observation = nobs, 
                    data = train_x)
explain_2
vd_xgb <- variable_importance(explainer_xgb, type = "raw")
plot(vd_xgb)



#混淆矩阵
pre_xgb = round(predict(model, newdata = xgb_test))
tab = table(test_y,pre_xgb,dnn=c("true","pre"))
#预测的准确率
caret::confusionMatrix(tab)


#ROC曲线
xgboost_roc <- roc(test_y,as.numeric(pre_xgb))
xgboost_roc
pdf('ROC.pdf',width =5 ,height =5 )
plot(xgboost_roc, print.auc=TRUE, auc.polygon=TRUE,
     grid=c(0.1, 0.2),grid.col=c("green", "red"),
     max.auc.polygon=TRUE,auc.polygon.col="skyblue",
     print.thres=TRUE,main='ROC curve')
dev.off()

#计算模型的误差
mean((test_y - pre_xgb)^2) #mse
caret::MAE(test_y, pre_xgb) #mae
caret::RMSE(test_y, pre_xgb) #rmse
