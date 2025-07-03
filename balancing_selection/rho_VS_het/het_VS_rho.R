setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/rho_VS_het")

# 第二步：读取平均杂合度文件
heterozygosity_file <- "./het/average_heterozygosity.txt"
heterozygosity_data <- read.table(heterozygosity_file, header = TRUE)

# 第三步：读取重组率文件并指定列名
recombination_file <- "./rho/rho.excludeROH.txt"
recombination_data <- read.table(recombination_file, header = FALSE, col.names = c("chromosome", "start", "end", "recombination_rate"))

# 为重组率数据添加杂合度列
recombination_data$average_heterozygosity <- NA

# 遍历重组率数据的每一行
for (i in 1:nrow(recombination_data)) {
  chr <- recombination_data$chromosome[i]
  start_rho <- recombination_data$start[i]
  end_rho <- recombination_data$end[i]
  
  # 筛选出当前重组率窗口所在的杂合度窗口
  relevant_het <- heterozygosity_data[heterozygosity_data$chromosome == chr & 
                                        heterozygosity_data$start <= start_rho & 
                                        heterozygosity_data$end >= end_rho, ]
  
  # 如果找到对应的杂合度窗口，取其杂合度值
  if (nrow(relevant_het) > 0) {
    recombination_data$average_heterozygosity[i] <- relevant_het$average_heterozygosity
  }
}

# 筛掉缺少杂合度或重组率的窗口
filtered_data <- recombination_data[!is.na(recombination_data$recombination_rate) & 
                                      !is.na(recombination_data$average_heterozygosity), ]

# 进行线性回归，此时 x 为杂合度，y 为重组率
model <- lm(recombination_rate ~ average_heterozygosity, data = filtered_data)
summary_model <- summary(model)
intercept <- coef(model)[1]
slope <- coef(model)[2]
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

# 打开 PDF 文件进行绘图保存
pdf("recombination_heterozygosity_plot.pdf", width = 8, height = 6)

# 绘制散点图和回归线
plot(filtered_data$average_heterozygosity, filtered_data$recombination_rate,
     main = "Recombination Rate vs Heterozygosity",
     xlab = "Average Heterozygosity",
     ylab = "Recombination Rate",
     pch = 16, col = "blue",
     cex = 0.6,
     cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.4)
abline(model, col = "red", lwd = 2)


# 在图中添加 P 值
legend("top", legend = paste("P =", format(p_value, digits = 3)), 
       bty = "n", cex = 1.0, text.col = "black", text.font = 1)

# 关闭 PDF 文件
dev.off()
