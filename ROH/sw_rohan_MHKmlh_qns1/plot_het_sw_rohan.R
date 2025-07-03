setwd("/data2/projects/zwang/m.hk/ROH/sw_rohan")

sw_het <- read.table(gzfile("mhk.50000.100.het.tsv.gz", 'rt'), header = F)
close(sw_het)
rohan_het <- read.table(gzfile('SZ-1.hEst.gz', 'rt'), header = F)
close(rohan_het)
sw_ROH <- read.table('mhk.hom.win.tsv', header = T)
rohan_ROH <- read.table('SZ-1.ROH.bed', header = F)

arrChr <- paste0('mhkscf_', c(1:23))



# 加载必要的绘图包
library(scales)
library(ggplot2)
library(gridExtra)

# 初始化存储绘图结果的列表
plots <- list()

# 循环处理每条染色体
for (chr in arrChr) {
#chr <- 'mhkscf_3'
  # 提取 sw 方法下当前染色体的 het 数据
  sw_het_chr <- sw_het[sw_het[, 1] == chr, ]
  sw_het_midpoint <- sw_het_chr[, 4]
  sw_het_heterozygosity <- sw_het_chr[, 9]
  
  # 提取 rohan 方法下当前染色体的 het 数据
  rohan_het_chr <- rohan_het[rohan_het[, 1] == chr, ]
  rohan_het_midpoint <- (rohan_het_chr[, 2] + rohan_het_chr[, 3]) / 2
  rohan_het_heterozygosity <- rohan_het_chr[, 5]
  
  # 提取 sw 方法下当前染色体的 ROH 数据
  sw_ROH_chr <- sw_ROH[sw_ROH[, 1] == chr, ]
  sw_ROH_start <- sw_ROH_chr[, 3]
  sw_ROH_end <- sw_ROH_chr[, 4]
  
  # 提取 rohan 方法下当前染色体的 ROH 数据
  rohan_ROH_chr <- rohan_ROH[rohan_ROH[, 1] == chr, ]
  rohan_ROH_start <- rohan_ROH_chr[, 2]
  rohan_ROH_end <- rohan_ROH_chr[, 3]
  
  # 创建 sw 方法绘图所需的数据框
  sw_df <- data.frame(
    Midpoint = sw_het_midpoint,
    Heterozygosity = sw_het_heterozygosity
  )
  
  # 创建 rohan 方法绘图所需的数据框
  rohan_df <- data.frame(
    Midpoint = rohan_het_midpoint,
    Heterozygosity = rohan_het_heterozygosity
  )
  
  # 找出 sw 方法杂合度的最大值
  sw_max_heterozygosity <- if (nrow(sw_df) > 0) max(sw_df$Heterozygosity) else 0
  # 计算 sw 方法水平横线的纵坐标
  sw_roh_y <- sw_max_heterozygosity / 10
  
  # 找出 rohan 方法杂合度的最大值
  rohan_max_heterozygosity <- if (nrow(rohan_df) > 0) max(rohan_df$Heterozygosity) else 0
  # 计算 rohan 方法水平横线的纵坐标
  rohan_roh_y <- rohan_max_heterozygosity / 10
  
  # 创建 sw 方法 ROH 绘图所需的数据框
  sw_ROH_df <- data.frame(
    xmin = sw_ROH_start,
    xmax = sw_ROH_end,
    y = -sw_roh_y
  )
  
  # 创建 rohan 方法 ROH 绘图所需的数据框
  rohan_ROH_df <- data.frame(
    xmin = rohan_ROH_start,
    xmax = rohan_ROH_end,
    y = -rohan_roh_y
  )
  
  # 绘制 sw 方法的曼哈顿图和 ROH 片段
  sw_plot <- ggplot() +
    geom_point(data = sw_df, aes(x = Midpoint, y = Heterozygosity), color = "blue", size = 0.2) +
    geom_segment(data = sw_ROH_df, aes(x = xmin, xend = xmax, y = y, yend = y), color = "black", size = 0.5) +
    labs(title = "SW", x = chr) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    )
#  sw_plot
  
  # 绘制 rohan 方法的曼哈顿图和 ROH 片段
  rohan_plot <- ggplot() +
    geom_point(data = rohan_df, aes(x = Midpoint, y = Heterozygosity), color = "green", size = 0.2) +
    geom_segment(data = rohan_ROH_df, aes(x = xmin, xend = xmax, y = y, yend = y), color = "black", size = 0.5) +
    labs(title = "Rohan", x = chr) +
    scale_y_continuous(labels = number_format(accuracy = 0.001)) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    )
  
  # 组合两个图
  combined_plot <- grid.arrange(sw_plot, rohan_plot, nrow = 2)
  
  # 将组合图添加到列表中
  plots[[chr]] <- combined_plot
}

# 确定网格布局的行数和列数
num_chromosomes <- length(arrChr)
num_cols <- min(num_chromosomes, 6)
num_rows <- ceiling(num_chromosomes / num_cols)

# 打开 PDF 文件
pdf("SZ-1.sw_rohan.het_ROH.pdf", width = 18, height = 12)

# 排列并打印所有图
do.call(grid.arrange, c(plots, nrow = num_rows, ncol = num_cols))

# 关闭 PDF 文件
dev.off()
