setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/rho_VS_het/het")


# 获取目录下所有后缀为 .bed 的文件
bed_files <- list.files(pattern = "\\.bed$", full.names = TRUE)

# 初始化一个空的数据框用于存储合并后的数据
combined_data <- data.frame()

# 循环读取每个 .bed 文件
for (file in bed_files) {
  # 读取文件
  data <- read.table(file, header = FALSE, col.names = c("chromosome", "start", "end", "heterozygosity"))
  
  # 将当前文件的数据添加到合并数据框中
  combined_data <- rbind(combined_data, data)
}

# 按染色体、起始位置和终止位置分组，计算杂合度的平均值
average_heterozygosity <- aggregate(heterozygosity ~ chromosome + start + end, data = combined_data, FUN = mean)

# 重命名列名
colnames(average_heterozygosity) <- c("chromosome", "start", "end", "average_heterozygosity")

# 按照染色体和窗口起始位置进行排序
average_heterozygosity <- average_heterozygosity[order(average_heterozygosity$chromosome, average_heterozygosity$start), ]

# 输出结果到文件
write.table(average_heterozygosity, file = "average_heterozygosity.txt", sep = "\t", na = "nan", quote = FALSE, row.names = FALSE)
