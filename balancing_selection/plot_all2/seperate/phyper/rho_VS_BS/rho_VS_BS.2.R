#-------------------------------------------------------------------------
setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/phyper/rho_VS_BS")

# 读取 10kb 窗口重组率文件
genome_10kb <- read.table("rho.excludeROH.txt", header = F, col.names = c("chromosome", "start", "end", "rho"))

# 读取三个种群组合的 50kb 倍数窗口文件
BS_files <- c("MHKjx.MHKmlh_qns.MHKjx.MHKmlh_qns.bothbalanced.bed", "MHKjx.MOPfq_pt.MHKjx.MOPfq_pt.bothbalanced.bed", "MHKmlh_qns.MOPfq_pt.MHKmlh_qns.MOPfq_pt.bothbalanced.bed")
BS_list <- lapply(BS_files, function(file) {
  read.table(file, header = F, col.names = c("chromosome", "start", "end"))
})

# 定义函数来计算 50kb 倍数窗口的重组率
check_BS_recombination <- function(BS_data, genome_data) {
  # 初始化一个空向量来存储每个 50kb 倍数窗口的重组率
  BS_rhos <- numeric(nrow(BS_data))
  
  # 遍历每个 50kb 倍数窗口
  for (i in seq_len(nrow(BS_data))) {
    BS_chrom <- BS_data$chromosome[i]
    BS_start <- BS_data$start[i]
    BS_end <- BS_data$end[i]
    
    # 找到被当前 50kb 倍数窗口包含的 10kb 窗口
    included_windows <- genome_data[genome_data$chromosome == BS_chrom & 
                                      genome_data$start >= BS_start & 
                                      genome_data$end <= BS_end, ]
    
    # 如果有包含的 10kb 窗口，则计算平均重组率
    if (nrow(included_windows) > 0) {
      BS_rhos[i] <- mean(included_windows$rho)
    } else {
      # 如果没有包含的 10kb 窗口，将重组率设为 NA
      BS_rhos[i] <- NA
    }
  }
  
  # 将重组率添加到 50kb 倍数窗口数据中
  BS_data$rho <- BS_rhos
  
  # 去除没有包含 10kb 窗口的行
  BS_data <- BS_data[!is.na(BS_data$rho), ]
  
  return(BS_data)
}

# 对每个种群组合进行判断
results <- lapply(BS_list, check_BS_recombination, genome_data = genome_10kb)

# 合并所有种群组合的结果
combined_BS <- do.call(rbind, results)

# 标记 BS 区域
#combined_BS$source <- "BS"

# 从全基因组数据中排除 BS 区域对应的 10kb 窗口
non_BS_data <- genome_10kb
for (i in seq_len(nrow(combined_BS))) {
  BS_chrom <- combined_BS$chromosome[i]
  BS_start <- combined_BS$start[i]
  BS_end <- combined_BS$end[i]
  non_BS_data <- non_BS_data[!(non_BS_data$chromosome == BS_chrom & 
                                   non_BS_data$start >= BS_start & 
                                   non_BS_data$end <= BS_end), ]
}

# 合并 BS 和 非BS 数据，标记为 BS + Non-BS
combined_all <- rbind(combined_BS, non_BS_data)
combined_all$source <- "BS + Non-BS"

# 标记 BS 区域
combined_BS$source <- "BS"

# 合并 BS 和 BS + 非BS 数据用于比较
final_data <- rbind(combined_BS, combined_all)


# 进行显著性比较
test_result <- wilcox.test(rho ~ source, data = final_data)
cat("全基因组与组合文件重组率的 Wilcoxon 秩和检验结果:\n")
print(test_result)

p_value <- test_result$p.value

# 加载 ggplot2 包
library(ggplot2)

# 绘制箱线图
p <- ggplot(final_data, aes(x = source, y = rho)) +
  geom_boxplot() +
  labs(y = "Recombination Rate", title = "Comparison of Recombination Rate between BS and All Regions") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 24),  # 调整标题字体大小
    axis.title = element_text(size = 24),  # 调整坐标轴标题字体大小
    axis.text = element_text(size = 24),   # 调整坐标轴刻度字体大小
    legend.text = element_text(size = 24)  # 调整图例字体大小（此图无图例，但可按需扩展）
  )

# 获取 y 轴的最大值
y_max <- max(final_data$rho, na.rm = TRUE)

# 计算横线的 y 位置
line_y <- y_max * 1.1

# 添加横线和 p 值
p <- p + 
  geom_segment(aes(x = 1, xend = 2, y = line_y, yend = line_y)) +
  annotate("text", x = 1.5, y = line_y * 1.05, label = paste("p =", format(p_value, digits = 3)), size = 8)

# 显示图形
print(p)
# 保存为 PDF 文件
ggsave("rho.BS_VS_All.pdf", plot = p, width = 10, height = 8)


# 计算并打印 BS 组的 95% 置信区间、平均数和中位数
BS_data <- final_data[final_data$source == "BS", "rho"]
BS_t_test <- t.test(BS_data)
BS_mean <- mean(BS_data)
BS_median <- median(BS_data)
cat("\nBS 组:\n")
cat("95% 置信区间:", BS_t_test$conf.int, "\n") #95% 置信区间: 0.1100522 0.4617493 
cat("平均数:", BS_mean, "\n") #平均数: 0.2859007 
cat("中位数:", BS_median, "\n") #中位数: 0.1138325 

# 计算并打印 BS + Non-BS 组的 95% 置信区间、平均数和中位数
BS_NonBS_data <- final_data[final_data$source == "BS + Non-BS", "rho"]
BS_NonBS_t_test <- t.test(BS_NonBS_data)
BS_NonBS_mean <- mean(BS_NonBS_data)
BS_NonBS_median <- median(BS_NonBS_data)
cat("\nBS + Non-BS 组:\n")
cat("95% 置信区间:", BS_NonBS_t_test$conf.int, "\n") #95% 置信区间: 0.1154421 0.118401 
cat("平均数:", BS_NonBS_mean, "\n") #平均数: 0.1169215 
cat("中位数:", BS_NonBS_median, "\n") #中位数: 0.07903201 
