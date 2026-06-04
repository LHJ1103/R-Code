library(limma)
library(survival)
library(survminer)
library(beepr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(jskm)
library(data.table)
library(doParallel)
setDTthreads(threads = detectCores())
KMdatas=read.csv("MYDATA.csv",header = TRUE)
KMdatas$Group=ifelse(KMdatas$gene>=median(KMdatas$gene),"Hight","Low")

diff=survdiff(Surv(futime,fustat)~Group,data = KMdatas)
pValue= 1-pchisq(diff$chisq,df=1)

rt_cut=surv_cutpoint(KMdatas,time = "futime",event = "fustat",variables = "gene")
rt_cat=surv_categorize(rt_cut)
colnames(rt_cat)=c("futime","fustat","gene")
diff=survdiff(Surv(futime,fustat)~ gene,data = rt_cat)
pValue= 1-pchisq(diff$chisq,df=1)
HR=(diff$obs[1]/diff$exp[1])/(diff$obs[2]/diff$exp[2])
up95=exp(log(HR)+qnorm(0.975*sqrt(1/diff$exp[2]+1/diff$exp[1])))
low95=exp(log(HR)-qnorm(0.975*sqrt(1/diff$exp[2]+1/diff$exp[1])))
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
gene="CREB5"
fit=surv_fit(Surv(futime,fustat)~gene,data = rt_cat)
surplot=ggsurvplot(fit,
                   data = KMdatas,
                   pval = paste0(pValue,"\n","HR=",round(HR,3)),
                   pval.size= 6,
                   surv.median.line = "hv",
                   legend.title= gene,
                   legend.labs=c("High level","Low level"),
                   xlab="Time(years)",
                   ylab= "Overall survival",
                   break.time.by=1,
                   risk.table = T,
                   risk.table.title="",
                   risk.table.col="strata",
                   risk.tables.height = 25)
