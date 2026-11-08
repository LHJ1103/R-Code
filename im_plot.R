
if (!require("ggcorrplot", quietly = TRUE)) {
  install.packages("ggcorrplot")
}


if (!require("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}


library(devtools)



library(CIBERSORT)
library(ggcorrplot)
library(tidyr)


DEG_expr <- read.csv("DEG_expr.csv",row.names = 1)

data <- data(LM22)

result <- cibersort(sig_matrix = LM22, 
                    mixture_file = DEG_expr, 
                    perm = 1, 
                    QN = TRUE 
                    )
class(DEG_expr)

DEG_expr <- as.matrix(DEG_expr)



cor_results <- read.csv("results.csv",row.names = 1,header = T,sep = ",")

head(cor_results)
table(cor_results$Gene)
table(cor_results$im_cell)

ggplot(cor_results, aes(x = Gene, y = im_cell, fill = Cor)) +   
  
  geom_tile() +   
  
  scale_fill_gradient2(low = "#7AA7EF", mid = "white", high = "#E9384E") +
  
  labs(title = "Correlation between immune cells and genes", x = "Gene", y = "Immune Cell") +  
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  
  geom_text(aes(label = ifelse(p.value < 0.05, format(p.value, nsmall = 2), "")), size = 3, color = "white", hjust = 0.5) +   
  
  scale_y_discrete(limits = rev(levels(cor_results$im_cell))) 




ggplot(cor_results, aes(x = Gene, y = im_cell, fill = Cor)) +
  geom_tile() +
  scale_fill_gradient2(low = "#7AA7EF", mid = "white", high = "#E9384E") +
  labs(title = "Correlation between immune cells and genes", x = "Gene", y = "Immune Cell") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  
  geom_text(aes(label = ifelse(p.value < 0.001, "***",  
                               ifelse(p.value < 0.01, "**",
                                      ifelse(p.value < 0.05, "*","")))), size = 3, color = "white", hjust = 0.5) +
  scale_y_discrete(limits = rev(levels(cor_results$im_cell))) +
  
  theme(axis.text.x = element_text(size = 15, face = "bold"),  
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title=element_text(size=20,face="bold"),
        panel.background = element_rect(fill = "white", color = "grey"),
        strip.background = element_rect(fill = "lightgray", color = "grey"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey", size = 0.5),
        panel.margin = unit(5, "lines"))  

  
  

  