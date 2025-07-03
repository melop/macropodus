setwd("/data2/projects/zwang/m.hk/ROH/simulation2/plot1")
# 设置要遍历的目录
directory_path <- "/data2/projects/zwang/m.hk/ROH/simulation2/plot1"
pattern <- "real"
# 获取目录中的文件列表
file_list <- list.files(directory_path, pattern = pattern)

# 设置保存PDF文件的路径和文件名
pdf_file_path <- "/data2/projects/zwang/m.hk/ROH/simulation2/plot1/all.pdf"


# 加载所需的包
library(ggplot2)

pdf(pdf_file_path, width = 20, height = 12)

# 循环遍历文件列表并生成图形
for (file_path in file_list) {
split_element <- strsplit(file_path, "\\.")[[1]]
groupname <- split_element[1]
# 读取数据文件
#print(groupname)

file1 <- paste(groupname,"intersect_len.txt", sep = ".")
file2 <- paste(groupname,"real.intersect_len.txt", sep = ".")
data <- scan(file1)
data2 <- read.table(file2, header = FALSE)

# 创建数据框
df <- data.frame(value = data)

# 绘制频率分布密度图
density_plot <- ggplot(df, aes(x = value)) +
  geom_density(fill = "skyblue", color = "black", alpha = 0.5) +
  labs(x = "SimulationLength", y = "Density") +
  ggtitle(groupname) +
  theme_minimal()+
  theme(
    text = element_text(size = 30)  # 设置所有文字的大小
  )
# 计算置信区间
confidence_level <- 0.95
mean_value <- mean(data)
se <- sd(data) / sqrt(length(data))
margin_of_error <- qt((1 + confidence_level) / 2, df = length(data) - 1) * se
confidence_interval <- c(mean_value - margin_of_error, mean_value + margin_of_error)

# 在密度图上添加置信区间的线
density_plot_with_ci <- density_plot +
  geom_vline(xintercept = confidence_interval, color = "red", linetype = "dashed") +
  annotate("text", x = mean(confidence_interval), y = 8e-07, label = "95% CI", color = "red", size = 10)

# 在密度图上添加来自data2的唯一值的垂直线
density_plot_with_ci_and_data2 <- density_plot_with_ci +
  geom_vline(aes(xintercept = data2$V1), color = "darkgreen", linetype = "dashed") +
  annotate("text", x = data2$V1, y = 8e-07, label = "RealLength", color = "darkgreen", size = 10)

# 调整X轴范围，使垂直线稍微左移
p <- density_plot_with_ci_and_data2 + coord_cartesian(xlim = c(min(data), max(data2$V1) + 1e+06))
print(p)
}
dev.off()