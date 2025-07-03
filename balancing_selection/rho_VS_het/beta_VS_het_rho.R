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
  
  # 提取染色体名称并转换格式，去掉种群名
  chrom_name <- gsub(".*chrom_", "mhkscf_", sub("^MOPdsy\\.", "", sub("\\.betascores.txt", "", basename(file))))
  
  # 添加染色体列
  data$chromosome <- chrom_name
  
  # 筛掉beta值小于0的碱基
  data <- data[data$beta_value >= 0, ]
  
  # 调整列顺序
  data <- data[, c("chromosome", "position", "beta_value")]
  
  # 将当前文件的数据添加到合并数据框中
  combined_data <- rbind(combined_data, data)
}

# 第二步：读取平均杂合度文件
# 假设平均杂合度文件在当前工作目录，这里需要根据实际情况调整路径
heterozygosity_file <- "./het/average_heterozygosity.txt"
heterozygosity_data <- read.table(heterozygosity_file, header = TRUE)

# 第三步：读取重组率文件并指定列名
recombination_file <- "./rho/rho.excludeROH.txt"
recombination_data <- read.table(recombination_file, header = FALSE, col.names = c("chromosome", "start", "end", "recombination_rate"))

# 第四步：为每个碱基分配杂合度和重组率
combined_data$heterozygosity <- NA
combined_data$recombination_rate <- NA
for (i in seq_len(nrow(combined_data))) {
  chrom <- combined_data$chromosome[i]
  pos <- combined_data$position[i]
  
  # 找到该碱基所在的杂合度窗口
  hetero_window <- heterozygosity_data[heterozygosity_data$chromosome == chrom & 
                                         heterozygosity_data$start <= pos & 
                                         heterozygosity_data$end >= pos, ]
  if (nrow(hetero_window) > 0) {
    combined_data$heterozygosity[i] <- hetero_window$average_heterozygosity
  }
  
  # 找到该碱基所在的重组率窗口
  rec_window <- recombination_data[recombination_data$chromosome == chrom & 
                                     recombination_data$start <= pos & 
                                     recombination_data$end >= pos, ]
  if (nrow(rec_window) > 0) {
    combined_data$recombination_rate[i] <- rec_window$recombination_rate
  }
}

# 第五步：进行线性回归，跳过缺失数据的碱基
complete_data <- combined_data[complete.cases(combined_data$beta_value, combined_data$heterozygosity, combined_data$recombination_rate), ]
regression_result <- lm(beta_value ~ heterozygosity + recombination_rate, data = complete_data)

# 输出回归结果
summary(regression_result)

summary_result <- summary(regression_result)
# 保存 summary 结果到文件
capture.output(summary_result, file = "regression_summary.txt")
#-----------------------------------------------------------------------------------------------
# 提取p值
p_value <- summary(regression_result)$coefficients[, 4]

# 打开PDF文件
pdf("beta_VS_het_rho.linear_regression_plot.pdf")

# 绘制散点图，这里以杂合度为x轴，beta值为y轴，用颜色区分不同重组率
plot(complete_data$heterozygosity, complete_data$beta_value, 
     xlab = "Heterozygosity", 
     ylab = "Beta Value", 
     main = "Linear Regression of Beta Value vs Heterozygosity and Recombination Rate",
     col = rainbow(length(unique(complete_data$recombination_rate)))[as.factor(complete_data$recombination_rate)])

# 添加回归直线（这里只是示意，实际可能需要更复杂的方式来展示三维关系）
abline(regression_result, col = "red")

# 在图中添加p值
legend("topright", legend = paste("p-values:\nHeterozygosity =", format(p_value["heterozygosity"], scientific = TRUE, digits = 2),
                                  "\nRecombination Rate =", format(p_value["recombination_rate"], scientific = TRUE, digits = 2)),
       bty = "n", text.col = "blue")

# 关闭PDF文件
dev.off()

# 输出回归结果
summary(regression_result)

# 保存合并后的数据到文件
write.table(combined_data, file = "combined_beta_heterozygosity_recombination.txt", sep = "\t", na = "nan", quote = FALSE, row.names = FALSE)

# 第六步：进行杂合度和重组率之间的线性回归
regression_result_hetero_rec <- lm(heterozygosity ~ recombination_rate, data = complete_data)

regression_result_rec_beta <- lm(beta_value ~ recombination_rate, data = complete_data)
regression_result_hetero_beta <- lm(beta_value ~ heterozygosity, data = complete_data)

regression_result_rec_residualbeta <- lm(regression_result_hetero_beta$residuals ~ complete_data$recombination_rate)

# 绘制散点图（杂合度和重组率关系）
plot(complete_data$recombination_rate, complete_data$heterozygosity, 
     xlab = "Recombination Rate", 
     ylab = "Heterozygosity", 
     main = "Linear Regression of Heterozygosity vs Recombination Rate")

# 添加回归直线（杂合度和重组率回归）
abline(regression_result_hetero_rec, col = "red")
