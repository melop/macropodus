setwd("/data2/projects/zwang/m.hk/ROH/simulation2/plot1")
# 设置要遍历的目录
directory_path <- "/data2/projects/zwang/m.hk/ROH/simulation2/plot1"
pattern <- "real"

# 获取目录中的文件列表
file_list <- list.files(directory_path, pattern = pattern)

# 设置保存PDF文件的路径和文件名
pdf_file_path <- "/data2/projects/zwang/m.hk/ROH/simulation2/plot1/Pvalue_hist.pdf"

# 创建空的向量存储每个循环的p值
p_values <- c()

pdf(pdf_file_path, width = 20, height = 12)

# 循环遍历文件列表并计算p值
for (file_path in file_list) {
  split_element <- strsplit(file_path, "\\.")[[1]]
  groupname <- split_element[1]
  
  # 读取数据文件
  file1 <- paste(groupname, "intersect_len.txt", sep = ".")
  file2 <- paste(groupname, "real.intersect_len.txt", sep = ".")
  data <- scan(file1)
  data2 <- read.table(file2, header = FALSE)$V1
  
  # 计算p值，即data中x值大于data2值的个数在data中总个数的比例
  p_value <- sum(data > data2) / length(data)
  
  # 将p值添加到向量中
  p_values <- c(p_values, p_value)
}

# 绘制频率分布直方图
hist_plot <- ggplot() +
  geom_histogram(aes(x = p_values, y = ..count..), binwidth = 0.002, fill = "skyblue", color = "black", boundary = 0) +
  labs(x = "P Value", y = "Count") +
#  ggtitle("Histogram of P Values") +
  theme_minimal()+
  xlim(0, 0.05)+  # 手动指定X轴范围
  ylim(0,100)+
  theme(
    text = element_text(size = 40),  # 设置文本字体大小
    axis.title = element_text(size = 40),  # 设置轴标题字体大小
#    plot.title = element_text(size = 14, hjust = 0.5)  # 设置图标题字体大小
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(), # 移除次要网格线
    axis.line = element_line(color = "black"),  # 设置轴线颜色
    axis.ticks = element_line(color = "black"),  # 设置刻度线颜色
    axis.ticks.length = unit(0.2, "cm"),  # 设置刻度线长度
    axis.text = element_text(color = "black")  # 设置刻度标签颜色
  )
# 打印频率分布直方图
print(hist_plot)

dev.off()
