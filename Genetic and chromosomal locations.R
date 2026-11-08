library(RCircos)
library(magrittr)
library(tidyverse)
library(rtracklayer)

genes <- c('EDN1','MMP7')

gene_pos <- import("gencode.v39.annotation.gtf") %>% 
  as.data.frame %>% 

  filter(source == "HAVANA", type == "gene") %>%  
  
  dplyr::select(seqnames, start, end, gene_name) %>%   
  
  filter(gene_name %in% genes)
head(gene_pos)


gene_pos$gene_dot <- rnorm(nrow(gene_pos), 0, 2)

write.csv(gene_pos,"very_easy_input.csv", row.names = F, quote = F)


gene_pos <- read.csv("very_easy_input.csv", header = T)


(data("UCSC.HG38.Human.CytoBandIdeogram"))



pdf(file="circGene.pdf", height=5, width=5)


RCircos.Set.Core.Components(UCSC.HG38.Human.CytoBandIdeogram,
                            chr.exclude = NULL,
                            tracks.inside = 3, 
                            tracks.outside = 0)
RCircos.Set.Plot.Area()


RCircos.Chromosome.Ideogram.Plot()


params <- RCircos.Get.Plot.Parameters()
params$track.background <- "grey" 
RCircos.Reset.Plot.Parameters(params)

RCircos.Scatter.Plot(gene_pos,
                     data.col = 5, 
                     by.fold = 1, 
                     track.num = 1,
                     side = "in")


RCircos.Gene.Connector.Plot(genomic.data = gene_pos, 
                            track.num = 2, 
                            side = "in")


RCircos.Gene.Name.Plot(gene_pos, 
                       name.col = 4,
                       track.num = 3, 
                       side = "in")

dev.off()




