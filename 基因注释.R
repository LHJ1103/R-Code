gene<-read.csv("GSE37283.csv",row.names = 1)
zhu<-read.csv("GPL13158-5065.csv")

#如果有多个基因，只取前面第一个
k<-zhu$Gene.Symbol
k<-strsplit(k,"///")
o<-data.frame()
for(i in 1:54613){
  k1<-k[[i]][1]
  o<-rbind(o,k1)
}
zhu$gene<-o$X.DDR1.


gene$ID<-rownames(gene)
m<-merge(zhu,gene,"ID")
m<-m[,-c(1:2)]
k=aggregate(.~gene,mean,data=m)
rownames(k)<-k$gene
k<-k[,-1]
write.csv(k,"GSE27283(注释).csv")
