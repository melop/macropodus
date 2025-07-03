# 设置工作目录
setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/phyper")

# 加载必要的包
library(VennDiagram)
library(grid)

# 定义函数用于统计文件中的区间数量
count_intervals <- function(file_path) {
  lines <- readLines(file_path)
  count <- length(lines)
  return(count)
}

# 指定要分析的两个种群
pop1 <- "MHKmlh_qns"
pop2 <- "MOPfq_pt"

# 构建文件名
file_both_balanced <- paste(pop1, pop2, pop1, pop2, "bothbalanced.bed", sep = ".")
file_only_sp1_balanced <- paste(pop1, pop2, "only", pop1, "balanced.bed", sep = ".")
file_only_sp2_balanced <- paste(pop1, pop2, "only", pop2, "balanced.bed", sep = ".")
file_both_no_balanced <- paste(pop1, pop2, pop1, pop2, "bothno.balanced.bed", sep = ".")

# 统计区间数量
k <- count_intervals(file_both_balanced)
m <- count_intervals(file_both_balanced) + count_intervals(file_only_sp1_balanced)
N_total <- count_intervals(file_both_balanced) + count_intervals(file_only_sp1_balanced) + count_intervals(file_only_sp2_balanced) + count_intervals(file_both_no_balanced)
n <- count_intervals(file_both_balanced) + count_intervals(file_only_sp2_balanced)

# 进行phyper检验
p_value <- phyper(k, m, N_total - m, n, lower.tail = FALSE)
print(paste("For combination", pop1, "and", pop2, "p - value:", p_value))

# 计算韦恩图所需的值
only_pop1 <- count_intervals(file_only_sp1_balanced)
only_pop2 <- count_intervals(file_only_sp2_balanced)
both <- count_intervals(file_both_balanced)
neither <- count_intervals(file_both_no_balanced)

# 打开PDF文件用于保存图形
#pdf(file = "Venn_diagram_with_rectangle_labeled.pdf", width = 10, height = 10)

# 创建一个新的绘图页面
grid.newpage()

# 绘制背景矩形
grid.rect(x = 0.5, y = 0.5, width = 0.9, height = 0.9, gp = gpar(fill = "lightgray", col = "black"))

# 绘制韦恩图
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.8, height = 0.8))
venn.plot <- draw.pairwise.venn(
  area1 = only_pop1 + both,
  area2 = only_pop2 + both,
  cross.area = both,
  category = c("", ""),
  fill = c("blue", "red"),
  alpha = 0.5,
  cat.pos = c(0, 0),
  cat.dist = c(0.05, 0.05),
  cex = 2,
  cat.cex = 0,
  lty = "blank",  # 去除圆的边框
  ext.text = TRUE,  # 允许外部文本
  ext.pos = 0,  # 外部文本位置
  ext.dist = 0.05  # 外部文本距离
)

# 添加标签到不同区域
# 两个种群都有的区域
if (length(venn.plot$intersect) > 0) {
  grid.text("MHKmlh_qns&MOPfq_pt with BSS\n", 
            x = venn.plot$intersect[1], y = venn.plot$intersect[2], gp = gpar(fontsize = 20))
} else {
  # 手动估算交集区域的大致坐标
  grid.text("MHKmlh_qns&MOPfq_pt with BSS\n",
            x = 0.6, y = 0.3, gp = gpar(fontsize = 20))
}
# 只有MHKmlh_qns的区域
if (length(venn.plot$circles[[1]]) > 0) {
  grid.text("only MHKmlh_qns BSS\n", 
            x = venn.plot$circles[[1]][1, 1], y = venn.plot$circles[[1]][1, 2], gp = gpar(fontsize = 20))
} else {
  # 手动估算只有MHKmlh_qns区域的大致坐标
  grid.text("only MHKmlh_qns BSS\n", 
            x = 0.3, y = 0.5, gp = gpar(fontsize = 20))
}
# 只有MOPfq_pt的区域
if (length(venn.plot$circles[[2]]) > 0) {
  grid.text("only MOPfq_pt BSS\n", 
            x = venn.plot$circles[[2]][1, 1], y = venn.plot$circles[[2]][1, 2], gp = gpar(fontsize = 20))
} else {
  # 手动估算只有MOPfq_pt区域的大致坐标
  grid.text("only MOPfq_pt BSS\n", 
            x = 0.85, y = 0.5, gp = gpar(fontsize = 20))
}

# 两个种群都不存在选择信号的区域（矩形内韦恩图外）
# 大致估算位置，这里可以根据实际情况调整
grid.text(paste("MHKmlh_qns&MOPfq_pt without BSS\n", neither), 
          x = 0.5, y = 0.05, gp = gpar(fontsize = 20))

popViewport()

# 关闭PDF文件
dev.off()

