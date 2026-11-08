
if (!require("e1071", quietly = TRUE)) {
  install.packages("e1071")
}


if (!require("caret", quietly = TRUE)) {
  install.packages("caret")
}

if (!require("pROC", quietly = TRUE)) {
  install.packages("pROC")
}


if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}


library(e1071)
library(caret)
library(pROC)
library(ggplot2)



mydata <- read.csv("data.csv",row.names = 1)


mydata <- subset(mydata, select = -time)

mydata$status <- factor(ifelse(mydata$status == "Tumor",1,0)) 

mydata$gender <- ifelse(mydata$gender == "female", 0, 1)
mydata$stage <- ifelse(mydata$stage == "early", 0, 1)
mydata$type <- ifelse(mydata$type == "squamous carcinom", 0, 1)
mydata$smoker <- ifelse(mydata$smoker == "nonsmoker", 0, 1)
mydata$alcohol_history <- ifelse(mydata$alcohol_history == "No", 0, 1)


table(mydata$status)
X <- mydata[,2:17] 
Y <- as.numeric(as.factor(mydata$status))

control <- trainControl(method = "repeatedcv",   
                        
                        number = 5,     
                        
                        repeats = 5,     
                        
                        search = "random"   
                        
                        )
  


set.seed(12345)   
svm_rfe <- rfe(X,
               Y,
               sizes = 1:16,
               rfeControl = rfeControl(functions = caretFuncs,
                                       
                                       method = "repeatedcv",
                                       
                                       number = 10,
                                       
                                       repeats = 5,
                                       
                                       verbose = FALSE),
                                      
               method = "svmLinear",
              
               trControl = control,
              
               preProc = c("center", "scale")
              
               )


svm_rfe_ranking <- svm_rfe$variables
head(svm_rfe_ranking)




varImp(svm_rfe)

varImp_dataframe <- data.frame(Gene = row.names(varImp(svm_rfe))[1:16],
                               importance = varImp(svm_rfe)[1:16, 1])
                          


varImp_dataframe <- na.omit(varImp_dataframe)


mycolors <- c('#D4E2A7','#88D7A4','#A136A1','#BAE8BC','#C757AF',
              '#DF9FCE','#D5E1F1','#305691','#B6C2E7','#E8EFF7',
              '#9FDFDF','#EEE0F5','#267336','#98CEDD','#CDE2EE',
              '#DAD490','#372E8A','#4C862D','#81D5B0','#BAE8C9',
              '#A7DCE2','#AFDE9C')

ggplot(varImp_dataframe, aes(x = reorder(Gene, -importance), y = importance , fill = Gene)) + 
   
  geom_col() +
  
  ggtitle("Hub Genes") +
  
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.title = element_text(margin = margin(b = 20)),
        panel.grid.major = element_line(color = "grey", size = 0.2)) +

  xlab("Gene") + ylab("Importance") +
  scale_fill_manual(values = mycolors)



top_10_vars <- svm_rfe_ranking$var[1:10]
top_10_vars


top10_SVM_data <- mydata[,top_10_vars]


X_plot = svm_rfe$results$Variables
Y_plot = svm_rfe$results$RMSE
plot(X_plot, Y_plot,  
     xlab="Variable Number",  
     ylab="RMSE (Cross-Validation)",  
     col="#7DEE44",   
     pch=16,                 
     cex=1.5,          
     lwd=2,            
     type="b",         
     ylim=c(0, 0.6)) 


abline(h=min(Y_plot), col="#DF294C")  
grid(col="grey",lwd=1,lty=3)   
 
legend(x = 6.8, y = 0.62,c("Training RMSE","Cross-Validation RMSE"),
       col=c("#7DEE44","#DF294C"),pch=c(16,NA),lwd=3,bg="white",cex = 1)   


wmin <- which.min(svm_rfe$results$RMSE)
wmin


points(wmin, svm_rfe$results$RMSE[wmin], col = "orange", pch = 16, cex=2)  
text(wmin, svm_rfe$results$RMSE[wmin],  
     paste0("N=", wmin), pos = 2, col = "orange", cex=1.5)  


Target_Genes <- svm_rfe$optVariables
Target_Genes


Best_SVM_data <- mydata[,Target_Genes]
write.csv(Best_SVM_data,"Best-svm-data.csv")


data <- cbind(mydata$status, Best_SVM_data)
colnames(data)[1] <- "status"

set.seed(12345)
train_index <- sample(1:nrow(data), nrow(data) * 0.7)
train_data <- data[train_index, ]
test_data <- data[-train_index, ]


model_linear <- svm(status~.,data = train_data,kernel="linear",probability = TRUE)


mean(train_data[,1] == model_linear$fitted)


table(actual = train_data[,1],model_linear$fitted) 


prob_linear <- predict(model_linear, train_data, probability=TRUE)
prob_linear_use <- attr(prob_linear,"probabilities")
linear_ROC <- roc(response = train_data$status,predictor = prob_linear_use[,2])


plot(linear_ROC,
     legacy.axes = TRUE,
     main="ROC curve",
     type = "l",col="red",lty=1,
     print.auc = T,
     thresholds="best", 
     print.thres="best") 





summary(data)
model_linear <- svm(status~.,data = train_data,kernel="linear",probability = TRUE)
model_polynomial <- svm(status~.,data = train_data,kernel="polynomial",probability = TRUE)
model_radial <- svm(status~.,data = train_data,kernel="radial",probability = TRUE)
model_sigmoid <- svm(status~.,data = train_data,kernel="sigmoid",probability = TRUE)

mean(train_data[,1] == model_linear$fitted)
mean(train_data[,1] == model_polynomial$fitted)
mean(train_data[,1] == model_radial$fitted)
mean(train_data[,1] == model_sigmoid$fitted)

table(actual = train_data[,1],model_linear$fitted) 
table(actual = train_data[,1],model_polynomial$fitted)
table(actual = train_data[,1],model_radial$fitted)
table(actual = train_data[,1],model_sigmoid$fitted)


accuracy <- c(mean(train_data[,1] == model_linear$fitted), 
              mean(train_data[,1] == model_polynomial$fitted),
              mean(train_data[,1] == model_radial$fitted),
              mean(train_data[,1] == model_sigmoid$fitted))
mycolors <- c('#54D9ED','#324B4F','#95B0B5','#D1BCFE','#FCFCD4','#F78ABA',
              '#00A9CB','#8EF8B6','#ED546C','#1BBEB9','#DE90DE','#00A44F')
barplot(accuracy, names.arg=c("Linear", "Polynomial", "Radial", "Sigmoid"), 
        main = "Accuracy of model", xlab = "Kernel Type", ylab = "Accuracy",
        col = mycolors)



prob_linear <- predict(model_linear, train_data, probability=TRUE)  
prob_linear_use <- attr(prob_linear,"probabilities")

prob_polynomial <- predict(model_polynomial, train_data, probability=TRUE)  
prob_polynomial_use <- attr(prob_polynomial,"probabilities")

prob_radial <- predict(model_radial, train_data, probability=TRUE)
prob_radial_use <- attr(prob_radial,"probabilities")  

prob_sigmoid <- predict(model_sigmoid, train_data, probability=TRUE)  
prob_sigmoid_use <- attr(prob_sigmoid,"probabilities")


linear_ROC <- roc(response = train_data$status, predictor = prob_linear_use[,2])  
polynomial_ROC <- roc(response = train_data$status, predictor = prob_polynomial_use[,2])
radial_ROC <- roc(response = train_data$status, predictor = prob_radial_use[,2])   
sigmoid_ROC <- roc(response = train_data$status, predictor = prob_sigmoid_use[,2])


linear_TPR <- linear_ROC$sensitivities  
linear_FPR <- 1-linear_ROC$specificities

polynomial_TPR <- polynomial_ROC$sensitivities 
polynomial_FPR <- 1-polynomial_ROC$specificities

radial_TPR <- radial_ROC$sensitivities
radial_FPR <- 1-radial_ROC$specificities

sigmoid_TPR <- sigmoid_ROC$sensitivities  
sigmoid_FPR <- 1-sigmoid_ROC$specificities  


ROC_df <- data.frame(FPR = c(linear_FPR, polynomial_FPR, radial_FPR, sigmoid_FPR), 
                     TPR = c(linear_TPR, polynomial_TPR, radial_TPR, sigmoid_TPR),
                     model = c(rep("Linear", length(linear_TPR)), 
                               rep("Polynomial", length(polynomial_TPR)),
                               rep("Radial", length(radial_TPR)),
                               rep("Sigmoid", length(sigmoid_TPR))))

 
plot <- ggplot(ROC_df, aes(x = FPR, y = TPR, col = model)) + 
  geom_line(size = 0.8) +
  geom_line(color = "#BC1328", linetype = 1, size = 0.8, data = subset(ROC_df, model == "Linear")) +     # Linear模型,红色线  
  geom_line(color = "#0067B5", linetype = 1, size = 0.8, data = subset(ROC_df, model == "Polynomial")) + # Polynomial模型,蓝色线
  geom_line(color = "#09891D", linetype = 1, size = 0.8, data = subset(ROC_df, model == "Radial")) +   # Radial模型,绿色线   
  geom_line(color = "orange", linetype = 1, size = 0.8, data = subset(ROC_df, model == "Sigmoid")) + # Sigmoid模型,橙色线
 
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +  
  
  ggtitle("ROC Curves") +    
  labs(x = "1-specificity", y = "Sensitivity") +
  
  theme(axis.text = element_text(face = "bold", size = 14, color = "black"),    
        axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
        axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0))))+
  
  theme(legend.position="bottom", legend.box = "horizontal")



plot + theme_grey() 
plot + theme_bw()  
plot + theme_classic() 