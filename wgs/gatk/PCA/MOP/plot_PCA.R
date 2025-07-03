setwd("/data2/projects/zwang/m.hk/Structure/PCA/MOP")
# 加载必要的库
library(ggplot2)
library(dplyr)

# 读取数据文件，假设数据文件名为data.txt，如果有不同的文件名，请相应修改
data_file <- "MOP.csv"
data <- read.table(data_file, header = FALSE, sep = ",")

# 提取数据中的PC1和PC2列
pca_data <- data[, c("V1", "V3", "V4")]  # 假设PC1在第2列，PC2在第3列

# 创建一个包含20种区分度明显的颜色的向量
distinct_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#1a1a1a", "#d11c24", "#356abc", "#c79810", "#755079",
  "#9a9a9a", "#005c96", "#f4cc70", "#ce8e00", "#b84600"
)

# 绘制PCA散点图
p <- ggplot(data = pca_data, aes(x = V3, y = V4, color = V1)) +
  geom_point(size = 3) +
  scale_color_manual(values = distinct_colors) +  # 使用自定义颜色
  labs(x = "PC 1", y = "PC 2", title = "MOP") +
  theme_minimal() +
  theme(
    text = element_text(size = 30),     # 调整文本字体大小
    legend.text = element_text(size = 30),  # 调整图例字体大小，这里设置为10
    panel.grid = element_blank(),  # 隐藏网格线
    axis.line = element_line(size = 0.5, color = "black"),  # 显示坐标轴，可以调整颜色和线条大小
    axis.title = element_text(size = 30)  # 调整坐标轴标题字体大小
  ) +
  scale_y_continuous(labels = scales::number_format(scale = 1, accuracy = 0.01)) +  # 设置纵坐标轴标签格式
  scale_x_continuous(labels = scales::number_format(scale = 1, accuracy = 0.01)) +  # 设置纵坐标轴标签格式
  labs(color = "Population")  # 设置图例标题

# 保存为PDF文件
pdf("PCA_MOP_refMHK.pdf",  width = 12, height = 10)
print(p)
dev.off()

