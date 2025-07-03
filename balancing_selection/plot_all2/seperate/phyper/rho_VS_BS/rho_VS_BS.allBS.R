#-------------------------------------------------------------------------
setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/phyper/rho_VS_BS")

# 读取 10kb 窗口重组率文件
#genome_10kb <- read.table("rho.excludeROH.txt", header = F, col.names = c("chromosome", "start", "end", "rho"))
genome_10kb <- read.table("rho.txt", header = F, col.names = c("chromosome", "start", "end", "rho"))

# 读取三个种群组合的 50kb 倍数窗口文件
#roh_files <- c("../../MHKjx.balanced.bed", "../../MHKmlh_qns.balanced.bed", "../../MOPfq_pt.balanced.bed", "../../MHKhn.balanced.bed", "../../MOPdxs_sg.balanced.bed", "../../MOPdsy.balanced.bed")
#roh_files <- c("../../MHKhn.balanced.bed", "../../MOPdxs_sg.balanced.bed", "../../MOPdsy.balanced.bed")
roh_files <- c("./MHKjx.MHKmlh_qns.MHKjx.MHKmlh_qns.bothbalanced.bed", "./MHKjx.MOPfq_pt.MHKjx.MOPfq_pt.bothbalanced.bed", "./MHKmlh_qns.MOPfq_pt.MHKmlh_qns.MOPfq_pt.bothbalanced.bed")

roh_list <- lapply(roh_files, function(file) {
#  read.table(file, header = F, col.names = c("chromosome", "start", "end", "BS_number"))
  read.table(file, header = F, col.names = c("chromosome", "start", "end"))
})

# 定义函数来计算 50kb 倍数窗口的重组率
check_roh_recombination <- function(roh_data, genome_data) {
  # 初始化一个空向量来存储每个 50kb 倍数窗口的重组率
  roh_rhos <- numeric(nrow(roh_data))
  
  # 遍历每个 50kb 倍数窗口
  for (i in seq_len(nrow(roh_data))) {
    roh_chrom <- roh_data$chromosome[i]
    roh_start <- roh_data$start[i]
    roh_end <- roh_data$end[i]
    
    # 找到被当前 50kb 倍数窗口包含的 10kb 窗口
    included_windows <- genome_data[genome_data$chromosome == roh_chrom & 
                                      genome_data$start >= roh_start & 
                                      genome_data$end <= roh_end, ]
    
    # 如果有包含的 10kb 窗口，则计算平均重组率
    if (nrow(included_windows) > 0) {
      roh_rhos[i] <- mean(included_windows$rho)
    } else {
      # 如果没有包含的 10kb 窗口，将重组率设为 NA
      roh_rhos[i] <- NA
    }
  }
  
  # 将重组率添加到 50kb 倍数窗口数据中
  roh_data$rho <- roh_rhos
  
  # 去除没有包含 10kb 窗口的行
  roh_data <- roh_data[!is.na(roh_data$rho), ]
  
  return(roh_data)
}

# 对每个种群组合进行判断
results <- lapply(roh_list, check_roh_recombination, genome_data = genome_10kb)

# 合并所有种群组合的结果
combined_roh <- do.call(rbind, results)

# 添加全基因组数据
genome_10kb$source <- "Whole Genome"
combined_roh$source <- "Sheared PR with BS"
#combined_roh <- combined_roh[,c(1:3,5,6)]
final_data <- rbind(genome_10kb, combined_roh)

y_max <- 1.0
y_min <- 0
# 绘制箱线图
boxplot(rho ~ source, data = final_data, 
        ylab = "Recombination Rate",
        ylim = c(y_min,y_max))

# 进行显著性比较
test_result <- wilcox.test(rho ~ source, data = final_data)
p_value <- test_result$p.value
# 计算横线的 y 坐标，设置在箱线图顶部稍微上方一点
y_pos <- y_max * 0.9  # 调整到合适的比例位置

# 在箱线图顶部添加横线
segments(x0 = 1, x1 = 2, y0 = y_pos, y1 = y_pos, lwd = 2)

# 在横线上方标注 P 值
text(x = 1.5, y = y_pos + 0.02 * (y_max - y_min), labels = paste("P =", format(p_value, scientific = TRUE, digits = 2)))

# 计算两个数据集的中位数
whole_genome_median <- median(final_data$rho[final_data$source == "Whole Genome"], na.rm = TRUE)
combined_roh_median <- median(final_data$rho[final_data$source == "Sheared PR with BS"], na.rm = TRUE)

# 输出中位数
cat("全基因组数据集的中位数是:", whole_genome_median, "\n")
cat("组合文件数据集的中位数是:", combined_roh_median, "\n")

