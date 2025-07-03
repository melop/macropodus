setwd("/data2/projects/zwang/m.hk/ROH/ROHan/test_MHKmls/liner_DP_VS_FROH")

# 读取表格数据
table1 <- read.csv("Froh_allDP4.MHKmls_test.csv", header = F, col.names = c("pop","ID","FROH"))
table2 <- read.csv("MHKmls.seq_depth.csv", header = F, col.names = c("RID","DP","ID"))


# 合并两个表格，基于个体名进行匹配
merged_table <- merge(table1, table2, by = "ID", all.x = FALSE)

# 按种群进行分组
populations <- unique(merged_table$pop)

# 创建一个PDF文件来保存所有的图形
pdf("FROH_DP_regression.MHKgbc_MHKlw.pdf", width = 10, height = 6)

# 循环对每个种群进行线性回归分析并绘图
for (pop in populations) {
  subset_data <- merged_table[merged_table$pop == pop, ]
  
  # 进行线性回归
  model <- lm(FROH ~ DP, data = subset_data)
  
  # 提取回归结果
  summary_model <- summary(model)
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  r_squared <- summary_model$r.squared
  p_value <- summary_model$coefficients[2, 4]
  
  # 绘制散点图和回归线
  plot(subset_data$DP, subset_data$FROH, 
       main = paste("pop:", pop, " - FROH_VS_DP"),
       xlab = "DP", ylab = "FROH",
       pch = 16, col = "blue")
  abline(model, col = "red", lwd = 2)
  
  # 在图中添加回归方程和R²、P值
  eq <- substitute(italic(y) == a + b %.% italic(x), list(a = format(intercept, digits = 2), b = format(slope, digits = 2)))
  legend("topleft", legend = c(as.character(eq), paste("R² =", format(r_squared, digits = 2)), paste("P =", format(p_value, digits = 3))), 
         bty = "n", cex = 0.8, text.col = c("black", "black", "black"), text.font = c(3, 1, 1))
}

# 关闭PDF文件
dev.off()
