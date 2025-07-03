#-------------------------------------------------------------------------------
setwd("/data2/projects/zwang/m.hk/ROH/plot_het_diffDP_diffFROH")

# 读取 SW 方法的数据
#sw_het <- read.table(gzfile("mhk.50000.100.het.tsv.gz", 'rt'), header = F)
#close(sw_het)
#sw_ROH <- read.table('test.1e-3.9.tsv', header = T)

# 定义样本名称
samples <- c("JX-2", "SZ-2", "FJ-1", "DSY-6", "DSY-7")

# 定义字体大小
font_size <- 20  # 可根据需求调整字体大小

# 加载必要的绘图包
library(scales)
library(ggplot2)
library(gridExtra)

# 打开 PDF 文件
pdf("diffFROH_diffDP.het_ROH.3.pdf", width = 12, height = 20)

# 定义染色体名称
arrChr <- paste0('mhkscf_', c(1:23))

# 定义绘图函数
plot_sample <- function(sample, chr) {
    # 读取当前样本的 rohan 方法数据
    rohan_het <- read.table(paste0(sample, ".rohan.het.bed"), header = F, sep = "\t")
    rohan_ROH <- read.table(paste0(sample, ".ROH.bed"), header = F)
    
    # 提取 rohan 方法下当前染色体的 het 数据
    rohan_het_chr <- rohan_het[rohan_het[, 1] == chr, ]
    rohan_het_midpoint <- (rohan_het_chr[, 2] + rohan_het_chr[, 3]) / 2
    rohan_het_heterozygosity <- rohan_het_chr[, 4]
    
    # 提取 rohan 方法下当前染色体的 ROH 数据
    rohan_ROH_chr <- rohan_ROH[rohan_ROH[, 1] == chr, ]
    rohan_ROH_start <- rohan_ROH_chr[, 2]
    rohan_ROH_end <- rohan_ROH_chr[, 3]
    
    het_data <- data.frame(
      Midpoint = rohan_het_midpoint,
      Heterozygosity = rohan_het_heterozygosity
    )
    
    roh_data <- data.frame(
      xmin = rohan_ROH_start,
      xmax = rohan_ROH_end,
      y = -0.0005
    )
    
    fill_color <- "green"
  
  
  # 找出杂合度的最大值
  max_heterozygosity <- if (nrow(het_data) > 0) max(het_data$Heterozygosity) else 0
  
  # 绘制图
  plot <- ggplot() +
    geom_area(data = het_data, aes(x = Midpoint, y = Heterozygosity), fill = fill_color, alpha = 1) +
    geom_line(data = het_data, aes(x = Midpoint, y = Heterozygosity), color = fill_color, size = 0.2) +
    geom_segment(data = roh_data, aes(x = xmin, xend = xmax, y = y, yend = y), color = "black", size = 1) +
    labs(title = sample, x = chr, y = "Heterozygosity") +
    coord_cartesian(ylim = c(-0.001, 0.01)) +  # 指定 y 轴范围
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = font_size),
      axis.title.y = element_text(size = font_size),
      axis.text = element_text(size = font_size),
      plot.title = element_text(size = font_size + 2),
      panel.grid = element_blank()
    )
  
  return(plot)
}

# 循环处理每条染色体
for (chr in arrChr) {
  # 存储每个样本的图
  plots <- list()
  
  # 循环处理每个样本
  for (sample in samples) {
    plot <- plot_sample(sample, chr)
    plots[[sample]] <- plot
  }
  
  # 组合所有样本的图
  combined_plot <- grid.arrange(grobs = plots, nrow = length(samples))
  
  # 打印组合图到 PDF
  print(combined_plot)
}

# 关闭 PDF 文件
dev.off()
