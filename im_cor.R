
if (!require("ggcorrplot", quietly = TRUE)) {
  install.packages("ggcorrplot")
}


if (!require("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}


library(ggcorrplot)
library(tidyr)


data <- read.csv("data.csv",row.names = 1)

cibersort_result <- read.csv("cibersort_result.csv",row.names = 1)

data <- t(DEG_expr)
cibersort_result <- cibersort_result[,1:22]


if(max(data) < 50) {data <- 2^data}


data <- t(DEG_expr)

identical(rownames(data),rownames(cibersort_result))


data <- data[,1:3]
data <- as.data.frame(data)

shapiro.test(data$SHISA3)

ks.test(data$LCP1, pnorm, mean = mean(data$LCP1), sd = sd(data$LCP1))

# 计算相关性 -------------------------------------------------------------------
colnames(data)
colnames(cibersort_result)
target_Gene <- "CREB5"
expr <- as.numeric(data[,target_Gene])

cor.test(cibersort_result$T.cells.CD8, expr, use = "everything", method = "spearman")


calculate_correlation <- function(Gene_expr, cibersort_result) {
  cor_matrix <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:ncol(cibersort_result)) {
    result <- cor.test(Gene_expr, cibersort_result[, i], method = "spearman")
    new_row <- data.frame("Gene" = "Gene", "im_cell" = colnames(cibersort_result)[i], "Cor" = result$estimate, "p-value" = result$p.value)
    cor_matrix <- rbind(cor_matrix, new_row)
  }
  return(cor_matrix)
}

calculate_correlation <- function(Gene_expr, cibersort_result) {
  cor_matrix <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:ncol(cibersort_result)) {
    result <- cor.test(Gene_expr, cibersort_result[, i], method = "pearson")
    new_row <- data.frame("Gene" = "Gene", "im_cell" = colnames(cibersort_result)[i], "Cor" = result$estimate, "p-value" = result$p.value)
    cor_matrix <- rbind(cor_matrix, new_row)
  }
  return(cor_matrix)
}


data_to_calculate <- data[,1:3]

results <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)

for (i in 1:ncol(data_to_calculate)) {
  print(i)
  gene_expr <- data_to_calculate[, i]
  corr_result <- calculate_correlation(gene_expr, cibersort_result)
  
  results <- rbind(results, corr_result)
}

colnames(data_to_calculate)

results$Gene <- c(rep("CREB5", 22),
                  rep("FPR1", 22),
                  rep("LCP1", 22))
                  
write.csv(results,"results.csv")
results <- results[,c(1,2,3)]
colnames(results)

results_wide <- spread(results, key = im_cell, value = Cor)

rownames(results_wide) <- results_wide$Gene

results_wide$Gene <- NULL


ggcorrplot(t(results_wide), 
           
           show.legend = T, 
           
           colors = c("#2166AC", "white", "#B2182B"), 
           
           digits = 2, 
           
           lab = T)
           
  



corr_mat <- cor(data, cibersort_result, method = "spearman")

ggcorrplot(corr_mat, 
           
           show.legend = T, 
           
           colors = c("#2166AC", "white", "#B2182B"), 
           
           digits = 2, 
           
           lab = T) + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
