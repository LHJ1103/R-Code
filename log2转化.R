# 读取CSV文件（假设第一列是基因名，作为行名）
tpm_data <- read.csv("GSE3629(注释).csv", row.names = 1)

# 检查数据前几行
head(tpm_data)

# 进行log2(n+1)转换
log_transformed <- log2(tpm_data + 1)

# 查看转换后的数据
head(log_transformed)

# 保存结果到新的CSV文件（可选）
write.csv(log_transformed, "GSE3629（log2_转化）.csv", quote = FALSE)

