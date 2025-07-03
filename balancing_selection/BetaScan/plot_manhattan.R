# 设置工作目录
setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/MOPdsy/beta_input_for_chrom")

# 加载所需的R包
library(ggplot2)
library(gridExtra)
library(ggpubr)

# 设置PDF输出文件名
pdf_output_file <- "chroms_manhattan1.pdf"

# 获取包含所有数据文件的文件列表
data_files <- list.files(pattern = "*.betascores.txt", full.names = TRUE)

# 读取区间文件
intervals <- read.table("/data2/projects/zwang/m.hk/ROH/ROHan/runs_of_het/MOPfq_pt/MOPfq_pt.het.intersect.bed", header = F, col.names = c("Chromosome", "Start", "End", "Het"), sep = "\t")


# 创建一个空的列表，用于存储所有曼哈顿图
manhattan_plots <- list()

# 循环处理每个数据文件
for (data_file in data_files) {
  # 从数据文件读取数据
  data <- read.table(data_file, header = TRUE, col.names = c("Position", "BetaValue"), sep = "\t")
  
  # 提取染色体名称
  chromosome_name <- gsub("MOPdsy.chrom_|\\.betascores\\.txt", "", basename(data_file))
  chromosome_name <- paste0("mhkscf_", chromosome_name)
  
  # 创建曼哈顿图
  manhattan_plot <- ggplot(data, aes(x = Position, y = BetaValue)) +
    geom_point(size = 0.5) +
    labs(title = basename(data_file)) +
    ylim(-30, 30)  # 在这里设置y轴范围
  
  # 读取与当前染色体匹配的区间
  matching_intervals <- intervals[intervals$Chromosome == chromosome_name, ]
  
  if (nrow(matching_intervals) > 0) {
    for (i in 1:nrow(matching_intervals)) {
      manhattan_plot <- manhattan_plot + geom_segment(data = matching_intervals[i, ], aes(x = Start, xend = End, y = 10, yend = 10), color = "blue")
    }
  }
  
  # 将曼哈顿图添加到列表中
  manhattan_plots <- c(manhattan_plots, list(manhattan_plot))
}

# 将所有曼哈顿图放在同一个PDF文件中，每个图占一页
pdf(pdf_output_file, width = 18, height = 10)
for (i in 1:length(manhattan_plots)) {
  print(manhattan_plots[[i]])
}
dev.off()




