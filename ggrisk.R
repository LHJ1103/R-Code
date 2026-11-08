
if (!require("ggrisk", quietly = TRUE)) {
  
  install.packages("ggrisk")
}

if (!require("survival", quietly = TRUE)) {
  
  install.packages("survival")
}

library(ggrisk)
library(survival)


mydata <- read.csv("ggrisk.csv",row.names = 1,header = T,sep = ",")


summary(mydata$time)
table(mydata$status)


mydata$time <- mydata$time/365


multi_cox_genes <- colnames(mydata)[3:5]
formula <- as.formula(paste("Surv(time, status) ~", paste(multi_cox_genes, collapse = " + ")))
formula
multi_cox_model <- coxph(formula, data = mydata)

summary(multi_cox_model)


ggrisk(multi_cox_model)

two_scatter(multi_cox_model)


ggrisk(multi_cox_model, 
       heatmap.genes = multi_cox_genes, 
       cutoff.value='median', 
       cutoff.x = 145, 
       cutoff.y = -0.8, 
       code.0 = 'Still Alive', 
       code.1 = 'Already Dead', 
       code.highrisk = 'High Risk', 
       code.lowrisk = 'Low Risk', 
       title.A.ylab='Risk Score', 
       title.B.ylab='Survival Time (year)', 
       title.A.legend='Risk Group', 
       title.B.legend='Status',
       title.C.legend='Expression', 
       size.ABC=1.5, 
       size.ylab.title=14, 
       size.Atext=11, 
       size.Btext=11, 
       size.Ctext=11, 
       size.yticks=0.5, 
       size.yline=0.5, 
       size.points=2, 
       size.dashline=1, 
       size.cutoff=5, 
       size.legendtitle=13, 
       size.legendtext=12, 
       color.A=c(low='#0591F1',high='#E12543'), 
       color.B=c(code.0='#0591F1',code.1='#E12543'), 
       color.C=c(low='#0591F1',median='white',high='#E12543'), 
       vjust.A.ylab=1, 
       vjust.B.ylab=2, 
       family='sans', 
       expand.x = 3,
       relative_heights = c(0.1,0.1,0.01,0.15))



two_scatter(multi_cox_model, 
            cutoff.value = 'median', 
            cutoff.x = 142, 
            cutoff.y = -0.5, 
            code.0 = 'Still Alive',
            code.1 = 'Dead', 
            code.highrisk = 'High Group', 
            code.lowrisk = 'Low Group', 
            title.A.legend = 'Riskscore', 
            title.B.legend = 'Event Status', 
            title.A.ylab = 'Riskscore', 
            title.B.ylab = 'Survival Time (year)', 
            title.xlab = 'Rank', 
            vjust.A.ylab = 1, 
            vjust.B.ylab = 2, 
            size.AB = 2, 
            size.ylab.title = 16, 
            size.xlab.title = 16, 
            size.Atext = 12, 
            size.Btext = 12,
            size.xtext = 12, 
            size.xyticks = 0.5, 
            size.xyline = 0.5, 
            size.dashline = 1.5, 
            size.points = 2,
            size.cutoff = 5, 
            size.legendtitle = 14,
            size.legendtext = 13, 
            color.A = c(low='#4B70E4',high='#E12543'), 
            color.B = c(code.0='#4B70E4',code.1='#E12543'), 
            family = 'sans', 
            expand.x = 10) 
