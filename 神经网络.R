data = read.csv("data.csv",row.names = 1)
library(neuralnet)
ind = sample(2,nrow(data),replace = T,prob = c(1,0))
trainset = data[ind==1,]
testset = data[ind==2,]
trainset$Tumor = trainset$fen.lasso=="Tumor"
trainset$UC = trainset$fen.lasso=="UC"

attach(data)
network=neuralnet(Tumor+UC~CREB5+FPR1+LCP1,trainset,hidden = 1)
plot(network)
network=neuralnet(Tumor+UC~CREB5+FPR1+LCP1,trainset,hidden = c(1,1))
plot(network)
network=neuralnet(Tumor+UC~CREB5+FPR1+LCP1,trainset,hidden = c(1,1,1))
plot(network)

network
network$result.matrix
head(network$generalized.weights[[1]])
plot(network)
par(mfrow=c(2,2))
gwplot(network,selected.covariate = "CREB5")
gwplot(network,selected.covariate = "FPR1")
gwplot(network,selected.covariate = "LCP1")

net.predict = compute(network,testset[-1])$net.result
net.predict

net.prediction=c("Tumor","UC")[apply(net.predict,1,which.max)]
predict.table=table(testset$fen.lasso,net.prediction)
predict.table
library(caret)
confusionMatrix(predict.table)

