getwd()
setwd("E:/Rdata/WGCNA/")
getwd()
library(WGCNA)
options(stringAsFactors = FALSE)
femData = read.csv("GSE3629(分析用).csv",row.names = 1)

dim(femData)
names(femData)
datExpro = data.frame(t(femData))
colnames(datExpro) = rownames(femData)
rownames(datExpro) = colnames(femData)

datExpr1<-datExpro
m.vars=apply(datExpro,2,var)
expro.upper=datExpro[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4])]
datExpro1<-data.matrix(expro.upper)

gsg = goodSamplesGenes(datExpro1,verbose = 3)
gsg$allOK

sampleTree = hclust(dist(datExpro1),method = "average")
plot(sampleTree,main = "Sample clustering to detect outliers",sub = "", xlab = "", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
abline(h = 100,col = "red")
clust = cutreeStatic(sampleTree,cutHeight = 1,minSize = 10)

table(clust)
keepSamples = (clust==0)
datExpr = datExpro1[keepSamples,]
nGenes = ncol(datExpr)
nGenes
nSamples = nrow(datExpr)
nSamples
dim(datExpr)
datExpr = as.data.frame(datExpr)
sampleTree_final = hclust(dist(datExpr),method = "average")
plot(sampleTree_final,main = "Sample clustering to detect outliers final",sub = "", xlab = "", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

#热图
data.clin = read.csv("ClinicalTrait.CSV",header = T)
rownames(data.clin) = data.clin$Sample_geo_accession
data.clin = data.clin[,-1]
data.clin = as.data.frame(t(data.clin))
data.clin$UC = ifelse(grepl("UC",data.clin$Sample_characteristics_ch1),1,0)
data.clin$Tumor = ifelse(grepl("Tumor",data.clin$Sample_characteristics_ch1),1,0)
data.clin = data.clin[,-1]
nGenes =ncol(datExpr)
nSamples = nrow(datExpr)
traitColors = numbers2colors(data.clin,signed = FALSE)
plotDendroAndColors(sampleTree_final,traitColors,
                    groupLabels = names(data.clin),
                    main = "Sample dendrogram and trait heatmap")

#软阈值
powers = c(c(1:10),seq(from = 12,to=20,by=2))
sft = pickSoftThreshold(datExpr,powerVector = powers,verbose = 5)
sft
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.9,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
sft$powerEstimate
net = blockwiseModules(datExpr,power = 10,
                       TOMType = "unsigned",minModuleSize = 30,
                       reassignThreshold = 0,mergeCutHeight = 0.25,
                       numericLabels = TRUE,pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)
table(net$colors)
mergedColors = labels2colors(net$colors)
table(mergedColors)
plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,hang = 0.03,
                    addGuide = TRUE,guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",sep = "."),quote = F)
}





moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic




MEs0 = moduleEigengenes(datExpr,moduleColors)$eigengenes
MEs = orderMEs(MEs0)
head(MEs)
moduleTraitCor = cor(MEs,data.clin,use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor,nSamples)
textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                   signif(moduleTraitPvalue,1),")",sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6,8.5,3,3));
pdf("E:/新Rdata/WGCNA/")

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(data.clin),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationship"))
pdf("module_trait_heatmap.pdf", width = 10, height = 8)
dev.off()
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(data.clin), yLabels = names(MEs), cex.lab = 1,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))



data.clin = as.data.frame(data.clin)
modNames = substring(names(MEs),3)
geneModuleMembership = as.data.frame(cor(datExpr,MEs,use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))
names(geneModuleMembership) = paste("MM",modNames,sep = "")
names(MMPvalue) = paste("p.MM",modNames,sep = "")
geneTraitSignificance = as.data.frame(cor(datExpr,data.clin,use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) = paste("GS.",colnames(data.clin),sep = "")
names(GSPvalue) = paste("p.GS.",colnames(data.clin),sep = "")
module = "green"
column = match(module,modNames)
turquoise_moduleGenes = names(net$colors)[which(moduleColors==module)]
MM = abs(geneModuleMembership[turquoise_moduleGenes,column])
GS = abs(geneTraitSignificance[turquoise_moduleGenes,1])
par(mfrow = c(1,1))
verboseScatterplot(MM,GS,xlab = paste("Module Membership in",module,"module"),
                   ylab = "Gene significance for Tumor",
                   main = paste("Module membership vs.gene significance\n"),
                   cex.main = 1.2,cex.lab = 1.2,cex.axis = 1.2,col = module)
                     
