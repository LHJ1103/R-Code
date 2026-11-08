
if (!require("glmnet", quietly = TRUE)) {
  install.packages("glmnet")
}


library(glmnet)


mydata <- read.csv("data.csv",row.names = 1)


mydata <- subset(mydata, select = -time)

colnames(mydata[,1:17])

mydata$status <- ifelse(mydata$status == "Tumor",1,0) 

mydata$gender <- ifelse(mydata$gender == "female", 0, 1)
mydata$stage <- ifelse(mydata$stage == "early", 0, 1)
mydata$type <- ifelse(mydata$type == "squamous carcinom", 0, 1)
mydata$smoker <- ifelse(mydata$smoker == "nonsmoker", 0, 1)
mydata$alcohol_history <- ifelse(mydata$alcohol_history == "No", 0, 1)



y <- as.matrix(mydata[, 1])  
x <- as.matrix(mydata[, 2:17]) 



set.seed(12345)
lasso_model <- glmnet(x, 
                      y, 
                      family = "binomial",
                      alpha = 1) 
print(lasso_model) 


plot(lasso_model,
     xvar = "lambda",
     label = F)

coef_lasso <- coef(lasso_model,
                   s = 0.071930) 
coef_lasso



cv_model <- cv.glmnet(x, y, family = "binomial",alpha = 1,nfolds = 10)


plot(cv_model)


lambda_min <- cv_model$lambda.min
lambda_min
lambda_1se <- cv_model$lambda.1se
lambda_1se


coef_cv <- coef(lasso_model, s = lambda_min)
coef_cv 


exp(coef_cv)



coef_cv <- as.matrix(coef_cv)
coef_cv <- data.frame(coef_cv)

coef_cv$OR <- exp(coef_cv$s1)
nonzero_vars <- rownames(coef_cv[coef_cv$OR != 1, ])
nonzero_vars <- nonzero_vars [2:6]


lasso_data <- mydata[,nonzero_vars]
