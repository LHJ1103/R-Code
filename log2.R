
tpm_data <- read.csv("GSE3629.csv", row.names = 1)

head(tpm_data)


log_transformed <- log2(tpm_data + 1)


head(log_transformed)


write.csv(log_transformed, "GSE3629（log2）.csv", quote = FALSE)

