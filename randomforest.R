
if (!require("randomForest", quietly = TRUE)) {

  install.packages("randomForest")
}


library(randomForest)


mydata <- read.csv("data.csv",row.names = 1)


mydata <- subset(mydata, select = -time)

colnames(mydata[,1:17])

mydata$status <- ifelse(mydata$status == "Tumor",1,0) 
mydata$gender <- factor(ifelse(mydata$gender == "female", 0, 1))
mydata$stage <- factor(ifelse(mydata$stage == "early", 0, 1))
mydata$type <- factor(ifelse(mydata$type == "squamous carcinom", 0, 1))
mydata$smoker <- factor(ifelse(mydata$smoker == "nonsmoker", 0, 1))
mydata$alcohol_history <- factor(ifelse(mydata$alcohol_history == "No", 0, 1))


y <- mydata[, 1]
x <- mydata[, -1]



rf_model <- randomForest(x, 
                         y, 
                         ntree = 500, 
                         importance = TRUE)

varImpPlot(rf_model, 
           type=1, 
           main="Variable Importance",
           n.var = 16,
           scale = T,
           cex = 1.2
           )


varImpPlot(rf_model, 
           type=2, 
           main="Variable Importance",
           n.var = 16,
           scale = T,
           cex = 1.2
)


varImpPlot(rf_model,n.var = 16,cex = 1.2)


importance <- importance(rf_model)

importance <- as.data.frame(importance)
top10_IncMSE <- head(rownames(importance[order(-importance$`%IncMSE`),]), 10)
top10_IncMSE


top10_ncNodePurity <- head(rownames(importance[order(-importance$IncNodePurity),]), 10)
top10_ncNodePurity

top10_IncMSE_data <- mydata[,top10_IncMSE]
top10_ncNodePurity_data <- mydata[,top10_ncNodePurity]


write.csv(top10_IncMSE_data,"top10_IncMSE_data.csv")
write.csv(top10_ncNodePurity_data,"top10_ncNodePurity_data.csv")
