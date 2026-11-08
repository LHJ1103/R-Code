
if (!require("VIM", quietly = TRUE)) {
  
  install.packages("VIM")
}


if (!require("rms", quietly = TRUE)) {
  
  install.packages("rms")
}


if (!require("nomogramFormula", quietly = TRUE)) {
  
  install.packages("nomogramFormula")
}


if (!require("pROC", quietly = TRUE)) {
  
  install.packages("pROC")
}


if (!require("rmda", quietly = TRUE)) {
  
  install.packages("rmda")
}


library(VIM)
library(rms)
library(nomogramFormula)
library(pROC)
library(rmda)
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


mydata <- read.csv("data.csv",row.names = 1,sep = ",",header = T)

val.data = read.csv("GSE37283YZ.csv",row.names = 1,sep = ",",header = T)
val.data2 = read.csv("GSE206285YZ.csv",row.names = 1,sep = ",",header = T)
val.data3 = read.csv("GSE75214YZ.csv",row.names = 1,sep = ",",header = T)
val.data4 = read.csv("GSE59071YZ.csv",row.names = 1,sep = ",",header = T)
mydata_copy <- mydata

aggr(mydata, col=c('skyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(mydata), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Percentage"))

mydata <- na.omit(mydata)
summary(mydata)

table(mydata$status)
mydata$status <- factor(mydata$status,level=c(0,1),labels=c("UC","Tumor"))



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
        dencol = "skyblue",
        boxcol = "skyblue")
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
     
     ROC$sensitivities,type = "l", 
     
     col="red", 
     
     lty=1,
     
     xlab = "1-Specificity",
     
     ylab = "Sensitivity",
     
     lwd=3) 

abline(0,1) 

legend(0.3,0.1,
       c("AUC of GSE206285 : 0.9595"),
       
       lty = c(1),   
       
       lwd = c(3),   
       
       col = c("red"), 
       
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
     
     ROC4$sensitivities,type = "l", 
     
     col="red", 
     
     lty=1,
     
     xlab = "1-Specificity",
     
     ylab = "Sensitivity",
     
     lwd=3) 

abline(0,1) 

legend(0.3,0.1,
       c("AUC of GSE59071 : 0.9597"),
       
       lty = c(1),   
      
       lwd = c(3),   
      
       col = c("red"), 
       
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


colnames(mydata)
formula <- as.formula(status ~ CREB5 + FPR1 + LCP1 )
model <- lrm(formula, 
           data = mydata, 
           x=TRUE, 
           y=TRUE) 
model

OR <- exp(model$coefficients)
OR


Nomogram_1 <- nomogram(model,
                       fun = function(x)1/(1+exp(-x)),
                       lp=F,
                       fun.at = c(0.1,0.5,0.9),
                       funlabel = "Risk of Tumor")
plot(Nomogram_1)



plot(Nomogram_1,
    
     xfrac = .35,
     
     cex.var = 1,
     
     cex.axis = 1,
     
     tcl = -0.5,
     
     lmgp = 0.3,
     
     label.every = 1,
     
     col.grid = gray(c(0.8,0.95)))






options(option) 

results <- formula_rd(nomogram = Nomogram_1)

val.data$points<-points_cal(formula = results$formula,rd=val.data)

head(mydata) 



install.packages("brglm") 
library(brglm)


model_ROC <- brglm(formula, data = mydata, family = binomial())
summary(model_ROC)

model_ROC <- glm(formula,data = mydata,family = binomial())



mydata$predvalue <- predict(model_ROC,type="response")

ROC <- roc(mydata$status,mydata$predvalue)
auc(ROC)
ci(auc(ROC))

plot(1-ROC$specificities,  
     
     ROC$sensitivities,type = "l", 
     
     col="red", 
     
     lty=1,
     
     xlab = "1-Specificity",
     
     ylab = "Sensitivity",
     
     lwd=3) 
     
abline(0,1) 

legend(0.5,0.25,
       c("AUC of Nomogram:0.8808"),
       
       lty = c(1),   
       
       lwd = c(3),   
       
       col = c("red"), 
       
       bty = "0")    


colnames(mydata)
roc1 <- roc(mydata$status,mydata$predvalue)
roc2 <- roc(mydata$status,mydata$CREB5)
roc3 <- roc(mydata$status,mydata$FPR1)
roc4 <- roc(mydata$status,mydata$LCP1)

plot(roc1, 
     
     col="red",  
    
     print.auc=F, 
     
     print.thres=F,
     
     lty=1,  
    
     xlab = "Specificity",
     
     ylab = "Sensitivity", 
     
     lwd=3)  
     
plot(roc2, add=TRUE, col="#99DD98",lwd=2.5)
plot(roc3, add=TRUE, col="#8F94E0",lwd=2.5)
plot(roc4, add=TRUE, col="#6A287E",lwd=2.5)
plot(roc5, add=TRUE, COL="#50EBEC")
plot(roc6, add=TRUE, COL="#FFFF00")
plot(roc7, add=TRUE, COL="#493D26")
plot(roc8, add=TRUE, COL="#C5908E")
plot(roc9, add=TRUE, COL="#6A287E")
plot(roc10, add=TRUE, COL="#7D0552")

auc(roc1)
ci(auc(roc1))

auc(roc2)
ci(auc(roc2))
auc(roc3)
ci(auc(roc3))
auc(roc4)
ci(auc(roc4))


legend(0.72,0.28,
       legend=c("Model:0.8808(95% CI: 0.7711-0.9905)",
                "CREB5:0.7558(95% CI: 0.5978-0.9138)",
                "FPR1:0.8169(95% CI: 0.6881-0.9456)",
                "LCP1:0.8343(95% CI: 0.7084-0.9602)"),
       col=c("red","#99DD98","#8F94E0","#6A287E"),
       lty=1,
       lwd=2,
       cex = 1)


cal <- calibrate(model,method ="boot" ,B=1000)


par(cex.lab = 1.2, cex.axis = 1.2, lwd = 2)


plot(cal, 
     xlim = c(0, 1),  
     ylim = c(0, 1), 
     xlab = "Predicted Probability",
     ylab = "Observed Probability", 
     subtitles = FALSE)


abline(0, 1, col = "red", lty = 4)


legend(0.001,1.1, 
       legend = c("Calibration Curve", "Ideal Line"), 
       col = c("black", "red"), 
       lty = c(1, 2), 
       lwd = 2)




plot(cal, 
     
     xlim=c(0,1),  
     
     ylim=c(0,1), 
     
     xlab="Predicted Probability",
     
     ylab="Observed Probability", 
     
     subtitles=FALSE) 
     




model_DCA <- decision_curve(formula, 
                           
                            data = mydata_copy, 
                            
                            family = binomial(link = 'logit'), 
                             
                            thresholds = seq(0,1,by=0.01),  
                            
                            confidence.intervals = 0.95, 
                           
                            study.design = 'case-control', 
                             
                            population.prevalence = 0.3)  
                            

plot_decision_curve(model_DCA,  
                    
                    curve.names = c('Nomogram'), 
                    
                    xlim = c(0,1.0),   
                    
                    cost.benefit.axis = FALSE,  
                    
                    col = c('red'),   
                    
                    confidence.intervals = FALSE, 
                    
                    standardize = FALSE) 
                    

plot_clinical_impact(model_DCA, 
                     
                     population.size = 1000,  
                     
                     cost.benefit.axis = T,    
                      
                     n.cost.benefits = 8, 
                     
                     col = c('red','blue'), 
                     
                     confidence.intervals = F) 
                     