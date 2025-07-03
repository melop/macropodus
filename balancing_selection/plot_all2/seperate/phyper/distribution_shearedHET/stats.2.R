setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/phyper/distribution_shearedHET")
library(ggplot2)

# 定义文件列表
file_list <- c("MHKjx", "MHKmlh_qns", "MOPfq_pt")

# 开启 PDF 图形设备
pdf("2pops.shearedHET.distribution.2.pdf", width = 24, height = 6)
    
# 遍历所有可能的文件对组合
for (i in seq_along(file_list)) {
  for (j in i:length(file_list)) {
    file_name <- paste0("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/", file_list[i], ".", file_list[j], ".intersected.HET.bed")
    if (file.exists(file_name)) {
      file <- read.table(file_name, header = F, sep = "\t", col.names = c("chr", "start", "end"))
      file$len <- file$end - file$start
          
      # 计算合适的分组边界，确保从 0 开始
      max_len <- max(file$len)
      breaks <- seq(0, max_len, by = 50000)
          
      # 使用 ggplot2 绘制直方图
      p <- ggplot(file, aes(x = len)) +
      geom_histogram(breaks = breaks, fill = "gray", color = "black") +
      labs(title = paste(file_list[i], "-", file_list[j]), x = "Length", y = "Frequency") +
      # 设置 x 轴范围为 0 - 1000000
      scale_x_continuous(expand = c(0, 0), limits = c(0, 1000000)) +
      # 设置 y 轴范围为 0 - 800
      scale_y_continuous(expand = c(0, 0), limits = c(0, 800)) +
      theme_minimal() +
      theme(
        # 设置标题字体大小
        plot.title = element_text(size = 24), 
        # 设置坐标轴标题字体大小
        axis.title = element_text(size = 24), 
        # 设置坐标轴刻度标签字体大小和颜色
        axis.text = element_text(size = 24, color = "black"), 
        # 设置图例标题字体大小
        legend.title = element_text(size = 24), 
        # 设置图例文本字体大小
        legend.text = element_text(size = 24) 
      )
          
      # 打印图形到 PDF
      print(p)
    }
  }
}
    
# 关闭 PDF 图形设备
dev.off()
    