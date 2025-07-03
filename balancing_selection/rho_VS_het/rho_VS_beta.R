# 设置工作目录，替换为实际目录
setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/rho_VS_het")

# 第一步：合并所有染色体文件
# 获取所有符合格式的文件
file_pattern <- "MOPdsy\\.chrom_\\d+\\.betascores.txt"
chrom_files <- list.files(path = "../MOPdsy/beta_input_for_chrom/", pattern = file_pattern, full.names = TRUE)

# 初始化一个空的数据框用于存储合并后的数据
combined_data <- data.frame()

# 循环读取每个文件
for (file in chrom_files) {
  # 读取文件的第二行到最后一行
  data <- read.table(file, header = FALSE, skip = 1, col.names = c("position", "beta_value"))
  
  # 筛掉 beta 值小于 0 的碱基
  data <- data[data$beta_value >= 0, ]
  
  # 提取染色体名称并转换格式，去掉种群名
  chrom_name <- gsub(".*chrom_", "mhkscf_", sub("^MOPdsy\\.", "", sub("\\.betascores.txt", "", basename(file))))
  
  # 添加染色体列
  data$chromosome <- chrom_name
  
  # 调整列顺序
  data <- data[, c("chromosome", "position", "beta_value")]
  
  # 将当前文件的数据添加到合并数据框中
  combined_data <- rbind(combined_data, data)
}

# 第二步：读取重组率文件
# 假设重组率文件在当前工作目录，这里需要根据实际情况调整路径
recombination_file <- "./rho/rho.excludeROH.txt"
recombination_data <- read.table(recombination_file, header = FALSE, col.names = c("chromosome", "start", "end", "recombination_rate"))

# 第三步：为每个碱基分配重组率
combined_data$recombination_rate <- NA
for (i in seq_len(nrow(combined_data))) {
  chrom <- combined_data$chromosome[i]
  pos <- combined_data$position[i]
  # 找到该碱基所在的窗口
  window <- recombination_data[recombination_data$chromosome == chrom & 
                                 recombination_data$start <= pos & 
                                 recombination_data$end >= pos, ]
  if (nrow(window) > 0) {
    combined_data$recombination_rate[i] <- window$recombination_rate
  }
}

# 第四步：进行线性回归，跳过缺失数据的碱基
complete_data <- combined_data[complete.cases(combined_data$beta_value, combined_data$recombination_rate), ]


# 进行线性回归
regression_result <- lm(beta_value ~ recombination_rate, data = complete_data)

# 提取 p 值
p_value <- summary(regression_result)$coefficients[2, 4]

# 打开 PDF 文件
pdf("rho_VS_beta.linear_regression_plot.pdf")

# 绘制散点图
plot(complete_data$recombination_rate, complete_data$beta_value, 
     xlab = "Recombination Rate", 
     ylab = "Beta Value", 
     main = "Linear Regression of Beta Value vs Recombination Rate",
     cex.lab = 1.5,  # 调整坐标轴标签字体大小
     cex.main = 1.5, # 调整主标题字体大小
     cex.axis = 1.5  # 调整坐标轴刻度字体大小
     )

# 添加回归直线
abline(regression_result, col = "red")

# 在图中添加 p 值
legend("topright", legend = paste("p-value =", format(p_value, scientific = TRUE, digits = 2)), 
       bty = "n", text.col = "blue")

# 关闭 PDF 文件
dev.off()

# 输出回归结果
summary(regression_result)

# 保存合并后的数据到文件
write.table(combined_data, file = "combined_beta_recombination.txt", sep = "\t", na = "nan", quote = FALSE, row.names = FALSE)
