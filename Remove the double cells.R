library(DoubletFinder)
#devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(Seurat)
library(dplyr)
library(future)
#BiocManager::install("glmGamPoi")
library(glmGamPoi)
#plan("multisession", workers = 4) #规划多线程
#options(future.globals.maxSize = 2000 * 1024^2)


rm_doublet <- function(name=NULL,label=NULL,input=NULL,dim.usage=30,auto="false") {
  EC <- readRDS(paste0(label,".rds"))
  dir.create(paste0("DoulbletFinder_",label))
  setwd(paste0("DoulbletFinder_",label))
  sink(paste0(label,"_log_doublet.txt"),type = c("output", "message"),append = TRUE, split=T)
  gc()
  inpath <- list.files(path = input,pattern = name,full.names = T)
  if (auto=="true") {
    inpath <- paste0(inpath,"/","04.Matrix/")
  }
  print(paste0("数据读取中！"," 开始时间:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
  
  # if (gene.column == "1") {
  #   EC <- Read10X(inpath, gene.column = 1)
  # } else {
  #   EC <- Read10X(inpath)
  # }
  # print(paste0(label,"读取完成: ",ncol(EC)," 个细胞","; ",nrow(EC)," 个基因"))
  # EC <- CreateSeuratObject(EC, project = label, min.features = 200)
  print(paste0(label,"构建对象: ",ncol(EC)," 个细胞","; ",nrow(EC)," 个基因"))
  EC$log10FeaturePerUMI <- log10(EC$nFeature_RNA) / log10(EC$nCount_RNA)
  EC <- PercentageFeatureSet(EC, pattern = "^mt-",col.name = "percent_mito")
  EC <- PercentageFeatureSet(EC, pattern = "^Hb[abdgmqz]|^Hbe[^g]",col.name = "percent_hb")
  EC <- subset(x = EC,
               subset = nFeature_RNA >= 200 &
                 nCount_RNA >= 500 &
                 log10FeaturePerUMI >= 0.8 &
                 percent_mito <= 15 &
                 percent_hb <= 3
  ) 
  print(paste0(label,"通过过滤: ",ncol(EC)," 个细胞","; ",nrow(EC)," 个基因"))
  
  print(paste0("数据预处理中！"," 当前时间:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
  ##去双胞需要先聚类，均使用默认参数即可，可以用管道符%>%写得简洁些
  # EC <- NormalizeData(EC)
  # EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 3000)
  # EC <- ScaleData(EC)
  EC <- SCTransform(EC,vst.flavor = "v2", verbose = T, method = "glmGamPoi",vars.to.regress = "percent_mito")
  EC <- RunPCA(EC)
  ElbowPlot(EC, ndims = 50)
  ggplot2::ggsave(filename = paste0(label,"_ElbowPlot.pdf"), width = 9, height = 7)
  EC <- RunUMAP(EC, dims = 1:dim.usage)
  EC <- RunTSNE(EC, dims = 1:dim.usage)
  EC <- FindNeighbors(EC, dims = 1:dim.usage) %>% FindClusters(resolution = 0.3)
  print(paste0(label,"的SCT矩阵: ",ncol(EC)," 个细胞","; ",nrow(EC)," 个基因"))
  print(paste0("数据预处理完成！"," 当前时间:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
  sink()
  p1 = DimPlot(EC, reduction = "umap",label = T)
  ggplot2::ggsave(p1,filename = paste0(label,"_doublet_UMAP_res0.3.pdf"), width = 9, height = 7)
  p2 = DimPlot(EC, reduction = "tsne",label = T)
  ggplot2::ggsave(p2,filename = paste0(label,"_doublet_tSNE_res0.3.pdf"), width = 9, height = 7)
  
  ##DoubletFinder去双胞的标准流程，封装成一个函数
  Find_doublet <- function(data){
    sink(paste0(label,"_log_doublet.txt"),type = c("output", "message"),append = TRUE, split=T)
    #寻找最优pK值
    print(paste0("正在寻找最优pK值！"," 当前时间:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
    sweep.res.list <- paramSweep(data, PCs = 1:dim.usage,sct = TRUE) ##使用log标准化，sct参数设置为 sct = F（默认 ）,如使用SCT标准化方法，设置为T
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    sink()
    pdf(paste0(label,"_find.pK.pdf"),width=4,height=4)
    bcmvn <- find.pK(sweep.stats) #可以看到最佳参数的点
    dev.off()
    sink(paste0(label,"_log_doublet.txt"),type = c("output", "message"),append = TRUE, split=T)
    p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK)) ##提取最佳pk值
    # 计算双细胞比例
    DoubletRate = ncol(data)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
    print(paste0("正在估算双细胞比例！"," 当前时间:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
    nExp_poi <- round(DoubletRate*ncol(data))  #nExp_poi <- round(0.05*ncol(data))
    #估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
    homotypic.prop <- modelHomotypic(data$seurat_clusters) #最好提供celltype，而不是seurat_clusters。
    # 使用同源双细胞比例对计算的双细胞比例进行校正 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    print(paste0("估算双细胞比例：",DoubletRate,sep = " "))
    print(paste0("估算双细胞数量：",nExp_poi,sep = " "))
    print(paste0("估算同源双细胞比例：",homotypic.prop,sep = " "))
    print(paste0("估算异源双细胞数量：",nExp_poi.adj,sep = " "))
    # 使用确定好的参数鉴定doublets
    print(paste0("正在运行doubletFinder鉴定doublets！"," 当前时间:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
    data <-
      doubletFinder(
        data,
        PCs = 1:dim.usage,
        pN = 0.25,
        pK = p,
        nExp = nExp_poi.adj,
        reuse.pANN = NULL,
        sct = TRUE
      )
    print(paste0("鉴定完成！"," 当前时间:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
    print(paste0("meta.data如下：",sep = " "))
    print(paste0(colnames(data@meta.data)))
    colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
    sink()
    return(data)
  }
  ##调用上面写好的函数，返回的是一个Seurat对象，meta.data信息里会有双胞信息，需要自己手动删除
  EC<-Find_doublet(EC)
  sink(paste0(label,"_log_doublet.txt"),type = c("output", "message"),append = TRUE, split=T)
  p1 = DimPlot(EC, reduction = "umap", group.by = "doublet_info")
  ggplot2::ggsave(p1,filename = paste0(label,"_doublet_UMAP.pdf"), width = 9, height = 7)
  p2 = DimPlot(EC, reduction = "tsne", group.by = "doublet_info")
  ggplot2::ggsave(p2,filename = paste0(label,"_doublet_tSNE.pdf"), width = 9, height = 7)
  
  EC@meta.data$library = label #顺便打上这个样本的label
  
  ##生成的Seurat对象有个问题，会在meta.data里多了很多pANN_开头的列，需要手动删除
  c <- grep("pANN_",colnames(EC@meta.data))
  EC@meta.data <- EC@meta.data[,-c]
  ##输出此样本的细胞数
  print(paste0(label,"运行完成！"," 结束时间:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
  sink()
  return(EC)
}


home <- "E:/UC-CRC/GSE264408UC/UC1-4/"
for (i in c("UC1")) {
  setwd(home)
  assign(i,rm_doublet(i,i,paste0(home,i,'/')))
}



colnames(UC1@meta.data)
VlnPlot(UC1,group.by = "doublet_info",
        features = c("nCount_RNA","nFeature_RNA"),
        pt.size = 0,ncol = 2)




setwd(home)
saveRDS(UC1,"UC1（未删除双细胞）.rds")


UC1 = subset(UC1,doublet_info=="Singlet")
saveRDS(UC1,"UC1（删除双细胞）.rds")





