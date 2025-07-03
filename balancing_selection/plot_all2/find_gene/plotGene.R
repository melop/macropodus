setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/find_gene")

library(ggplot2)

# 读取数据
data <- read.table("balanced.gene.txt", stringsAsFactors = FALSE)
# 手动添加列名
colnames(data) <- c("Population", "GeneList")

# 创建一个空的数据框用于存储处理后的数据
result <- data.frame(Population = character(), Gene = character(), stringsAsFactors = FALSE)

#-----------------------------------------------------------------------------------------------------------------
# 循环处理每一行，将基因拆分并存储
#for (i in 1:nrow(data)) {
#  genes <- unlist(strsplit(data$GeneList[i], ","))  # 拆分基因
#  populations <- rep(data$Population[i], length(genes))  # 重复种群名
#  result <- rbind(result, data.frame(Population = populations, Gene = genes, stringsAsFactors = FALSE))  # 合并
#}
#----------------------------------------------------------------------------------------------------------------

# 循环处理每一行，将组合的种群名拆分，并将基因对应到每个单独的种群
for (i in 1:nrow(data)) {
	  populations <- unlist(strsplit(data$Population[i], "\\."))  # 拆分种群名
  genes <- unlist(strsplit(data$GeneList[i], ","))  # 拆分基因列表
    
    # 为每个单独的种群创建对应的基因数据
    for (pop in populations) {
	        result <- rbind(result, data.frame(Population = pop, Gene = genes, stringsAsFactors = FALSE))
    }
}

# 为每个基因赋值1，以便在气泡图中显示
result$value <- 1


# 绘制气泡图
my_plot <- ggplot(result, aes(x = Gene, y = Population)) +
  geom_point(aes(size = value), shape = 21, fill = "blue", color = "black") +
  scale_size_continuous(range = c(5, 10)) +  # 调整气泡大小范围
  labs(x = "Gene", y = "Population") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), # 旋转x轴标签以提高可读性
  plot.title = element_text(size = 24),        # 设置标题字体大小
      axis.title = element_text(size = 28),        # 设置坐标轴标题字体大小
      axis.text = element_text(size = 28, color = "black"),         # 设置坐标轴刻度字体大小
      legend.title = element_text(size = 20),      # 设置图例标题字体大小
      legend.text = element_text(size = 20))        # 设置图例文本字体大小)
print(my_plot)

# 使用 ggsave 保存图形
ggsave("plot_output.2.pdf", plot = my_plot, width = 18, height = 6)
