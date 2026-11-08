if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)


getGDCprojects()$project_id 
project="TCGA-COAD"
TCGA_COAD <- GDCquery(project = project,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(query = TCGA_COAD,method="api")

expr <- GDCprepare(query=TCGA_COAD)

names(expr@assays)


counts <- as.data.frame(assay(expr)) 
TPM <- as.data.frame(assay(expr,i = "tpm_unstrand")) 


data=as.data.frame(rowRanges(expr))
colnames(data)
mydata <- data[,c("gene_type","gene_name")]
table(mydata$gene_type)
expr_count <-  cbind(gene_type=data$gene_type,gene_name=data$gene_name,counts)
expr_TPM <- cbind(gene_type=data$gene_type,gene_name=data$gene_name,TPM)


project = "TCGA-COAD"
query <- GDCquery(project = project,
                  data.category = "Clinical",
                  data.type = "Clinical Supplement")
                


GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")


save(expr_count,expr_TPM,clinical,file = "TCGA_COAD.Rdata")


load("TCGA_COAD.Rdata")


library(dplyr) 
table(expr_count$gene_type)
mRNA <- expr_count[expr_count$gene_type=="protein_coding",]

mRNA <- mRNA[,-1]
dim(mRNA)
mRNA <- mRNA %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), mean))


mRNA$mean <- rowMeans(mRNA[,2:523])
mRNA <- arrange(mRNA,desc(mean))
mRNA <- mRNA %>% distinct(gene_name, .keep_all = T)
mRNA <- select(mRNA,-c(mean))



mRNA <- as.data.frame(mRNA)
rownames(mRNA) <- mRNA$gene_name
mRNA <- mRNA[,-1]

table(substr(colnames(mRNA),14,16))

Tumor <- grep('01A',colnames(mRNA))
Tumor 
Tumor_mRNA <- mRNA[,Tumor]

Normal <- grep('11A',colnames(mRNA))
Normal 
Normal_mRNA <- mRNA[,Normal]

mRNA_count <- cbind(Normal_mRNA,Tumor_mRNA)
write.csv(mRNA_count,file = "mRNA_count.csv")


mRNA_TPM <- expr_TPM[expr_TPM$gene_type=="protein_coding",]
mRNA_TPM <- mRNA_TPM[,-1]
dim(mRNA_TPM)
mRNA_TPM <- mRNA_TPM %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), mean))
mRNA_TPM <- as.data.frame(mRNA_TPM)



mRNA_TPM$mean <- rowMeans(mRNA_TPM[,2:175])
mRNA_TPM <- arrange(mRNA_TPM,desc(mean))
mRNA_TPM <- mRNA_TPM %>% distinct(gene_name, .keep_all = T)
mRNA_TPM <- select(mRNA_TPM,-c(mean))

rownames(mRNA_TPM) <- mRNA_TPM$gene_name
mRNA_TPM <- mRNA_TPM[,-1]
table(substr(colnames(mRNA_TPM),14,16))

Tumor <- grep('01A',colnames(mRNA_TPM))
Tumor 
Tumor_mRNA_TPM <- mRNA_TPM[,Tumor]

Normal <- grep('11A',colnames(mRNA))
Normal 
Normal_mRNA_TPM <- mRNA_TPM[,Normal]

mRNA_TPM_new <- cbind(Normal_mRNA_TPM,Tumor_mRNA_TPM)
write.csv(mRNA_TPM,"mRNA_TPM.csv")

save(mRNA_count,mRNA_TPM_new,file = "count_and_TPM.Rdata")
