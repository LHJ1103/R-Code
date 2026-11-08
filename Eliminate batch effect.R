library(sva)
library(limma)
mm<-read.csv("GSE37283(UC+UCD).csv",row.names = 1)
mm2<-read.csv("GSE3629（log2）.csv",row.names = 1)

GSE37283_exp<-mm
boxplot(GSE37283_exp[,1:15],las="2")
GSE37283_group<-data.frame(row.names = colnames(GSE37283_exp),
                       Sample=colnames(GSE37283_exp),
                       DataSet="GSE37283")
GSE3629_exp<-mm2
boxplot(GSE3629_exp[,1:59],las="2")
GSE3629_group<-data.frame(row.names = colnames(GSE3629_exp),
                       Sample=colnames(GSE3629_exp),
                       DataSet="GSE3629")

dat_group<-rbind(GSE37283_group,GSE3629_group)

com_ensg<-intersect(rownames(GSE37283_exp),rownames(GSE3629_exp))

dat_exp_before<-cbind(GSE37283_exp[com_ensg,],
                      GSE3629_exp[com_ensg,])

boxplot(dat_exp_before[,1:74],las="2")

library(factoextra)
library(FactoMineR)
before_exp_pca<-PCA(t(dat_exp_before[,rownames(dat_group)]),
                    scale.unit=T,ncp=5,graph=F)
before_exp_pca.plot<-fviz_pca_ind(before_exp_pca,
                                  axes=c(1,2),
                                  label="none",
                                  addEllipses = T,
                                  ellipse.level=0.9,
                                  habillage = factor(dat_group$DataSet),
                                  palette = "aaas",
                                  mean.point=F,
                                  title="")
before_exp_pca.plot
before_exp_pca.plot<-before_exp_pca.plot+theme_bw()+
  theme(legend.direction = "horizontal",legend.position = "top")+
  xlim(-150,150)+ylim(-150,150)+
  xlab("Dim1")+ylab("Dim2")
before_exp_pca.plot


library(limma)
dat_exp<-removeBatchEffect(dat_exp_before[,rownames(dat_group)],
                           batch = dat_group$DataSet)
dim(dat_exp)
boxplot(dat_exp[,1:74],las="2")
dat_exp<-normalizeBetweenArrays(dat_exp)
boxplot(dat_exp[,1:74],las="2")
write.csv(dat_exp,"dat_exp.csv")


afterexp_pca<-PCA(t(dat_exp[,rownames(dat_group)]),
                    scale.unit=T,ncp=5,graph=F)
afterexp_pca.plot<-fviz_pca_ind(afterexp_pca,
                                  axes=c(1,2),
                                  label="none",
                                  addEllipses = T,
                                  ellipse.level=0.9,
                                  habillage = factor(dat_group$DataSet),
                                  palette = "aaas",
                                  mean.point=F,
                                  title="")
afterexp_pca.plot
afterexp_pca.plot<-afterexp_pca.plot+theme_bw()+
  theme(legend.direction = "horizontal",legend.position = "top")+
  xlim(-150,150)+ylim(-150,150)+
  xlab("Dim1")+ylab("Dim2")
afterexp_pca.plot

dat_PCA<-cowplot::plot_grid(before_exp_pca.plot,
                            afterexp_pca.plot,
                            ncol=2,nrow=1,
                            labels=toupper(letters)[1:2],
                            align = "hv")
dat_PCA
ggsave(plot=dat_PCA,
       "dat_PCA.pdf",
       width = 10,heigh=5,device=cairo_pdf)

getwd()
