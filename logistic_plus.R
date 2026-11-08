
if (!require("rms", quietly = TRUE)) {
  install.packages("rms")
}


if (!require("regplot", quietly = TRUE)) {
  install.packages("regplot")
}


if (!require("DynNom", quietly = TRUE)) {
  install.packages("DynNom")
}


library(rms) 
library(regplot) 
library(DynNom) 


mydata <- read.csv("data.csv",row.names = 1,sep = ",",header = T)

mydata$status <- factor(mydata$status,level=c(0,1),labels=c("UC","Tumor"))
mydata$gender <- factor(mydata$gender,level=c(0,1),labels=c("female","male"))
mydata$stage <- factor(mydata$stage,level=c(0,1),labels=c("early","advanced"))


dd = datadist(mydata)
option <- options(datadist = "dd")


colnames(mydata)
formula <- as.formula(status ~ CREB5 + FPR1 + LCP1 )
model <- lrm(formula,
             data = mydata, 
             x=TRUE, 
             y=TRUE) 


Nomogram_1 <- nomogram(model,
                       fun = function(x)1/(1+exp(-x)), 
                       lp=F, 
                       fun.at = c(0.1,0.5,0.9), 
                       funlabel = "Risk of Tumor")

plot(Nomogram_1,
     xfrac = .30,
     cex.var = 1.6,
     cex.axis = 1.4,
     tcl = -0.5,
     lmgp = 0.3,
     label.every = 1,
     col.grid = c("lightblue","#7CEF7A")
     )

Nomogram_2 <- glm(formula, family=binomial(link = logit), data=mydata)


regplot(Nomogram_2, 
        plots = c("density","boxes") ,
        observation = mydata[5,],
        center = T, 
        points = TRUE,
        droplines = F,
        title = "individual nomogram",
        odds = F, 
        interval="confidence",
        rank = NULL,
        clickable = F,
        dencol = "skyblue",
        boxcol = "skyblue"
        )



Nomogram_3 <- glm(formula,data = mydata, family = binomial)

DynNom(Nomogram_3,
       clevel = 0.55,
       DNtitle = "Nomogram",
       DNxlab = "probability",
       DNylab = NULL, 
       DNlimits = NULL
       )
