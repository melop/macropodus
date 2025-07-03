setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/phyper/distance_between_BS")
dat <- read.table("gff.genes.length.txt", header = F, sep = "\t")

MHKjx_MHKmlh_qns <- read.table("MHKjx.MHKmlh_qns.BS_on_shearedHET.minDistance.txt", header = T, sep = "\t")
MHKjx_MOPfq_pt <- read.table("MHKjx.MOPfq_pt.BS_on_shearedHET.minDistance.txt", header = T, sep = "\t")
MHKmlh_qns_MOPfq_pt <- read.table("MHKmlh_qns.MOPfq_pt.BS_on_shearedHET.minDistance.txt", header = T, sep = "\t")

# 定义字体大小
font_size <- 24  # 可根据需求调整字体大小

# 收集所有用于确定 x 轴范围的数据
all_x_values <- c(dat$V4, MHKjx_MHKmlh_qns$minDistance, MHKjx_MOPfq_pt$minDistance, MHKmlh_qns_MOPfq_pt$minDistance)

dat_x_max <- max(dat$V4)
breaks <- seq(0, dat_x_max, by = 10000)

# 计算 x 轴的最小值和最大值
x_min <- min(all_x_values)
x_max <- max(all_x_values)

# 计算 dat 第四列数据的 99% 置信区间
conf_int <- quantile(dat$V4, probs = c(0.005, 0.995))

#计算基因平均长度
filtered_values <- dat[, 4][dat[, 4] > 302 & dat[, 4] < 114959]
# 计算筛选后值的平均值
average_value <- mean(filtered_values)

# 加载 ggplot2 包
library(ggplot2)

# 开启 PDF 图形设备，设置宽度和高度（单位：英寸）
pdf("BSDistance_vs_geneLength.with99CI.pdf", width = 24, height = 8)

# 绘制直方图，设置组距为 10000
p <- ggplot(dat, aes(x = V4)) +
  geom_histogram(breaks = breaks, fill = "gray", color = "black") +
  labs(x = "Gene Length", y = "Frequency", title = "Histogram with Vertical Lines") +
  coord_cartesian(xlim = c(0, x_max)) +  # 手动设置 x 轴范围
  theme_minimal() +
  theme(
    axis.title = element_text(size = font_size),  # 坐标轴标题字体大小
    axis.text = element_text(size = font_size, color = "black"),   # 坐标轴刻度字体大小
    plot.title = element_text(size = font_size + 2),  # 图形标题字体大小，稍大一点
    legend.title = element_text(size = font_size),  # 图例标题字体大小
    legend.text = element_text(size = font_size)    # 图例文本字体大小
  )

# 定义绘制垂直线的函数
draw_vertical_lines <- function(data, color) {
  lapply(data$minDistance, function(x) {
    annotation_custom(
      grid::linesGrob(gp = grid::gpar(col = color, lty = 2, lwd = 2)),
      xmin = x, xmax = x, ymin = -Inf, ymax = Inf
    )
  })
}

# 添加垂直线
p <- p + 
  draw_vertical_lines(MHKjx_MHKmlh_qns, "red") +
  draw_vertical_lines(MHKjx_MOPfq_pt, "blue") +
  draw_vertical_lines(MHKmlh_qns_MOPfq_pt, "darkgreen")

# 添加 95% 置信区间的垂直线
p <- p + 
  geom_vline(xintercept = conf_int[1], color = "black", linetype = "dashed", size = 1) +
  geom_vline(xintercept = conf_int[2], color = "black", linetype = "dashed", size = 1)

# 创建虚拟数据用于生成图例
legend_data <- data.frame(
  label = c("MHKjx - MHKmlh_qns", "MHKjx - MOPfq_pt", "MHKmlh_qns - MOPfq_pt", "95% CI"),
  color = c("red", "blue", "darkgreen", "black")
)

# 添加虚拟的几何对象来触发图例生成
p <- p + geom_segment(
  data = legend_data,
  aes(x = -Inf, xend = -Inf, y = -Inf, yend = -Inf, color = label),
  linetype = 2,  # 设置虚线类型
  size = 2,
  show.legend = TRUE
) +
  scale_color_manual(
    name = "Sample Pairs",
    values = setNames(legend_data$color, legend_data$label)
  )

# 显示图形
print(p)
# 关闭 PDF 图形设备
dev.off()