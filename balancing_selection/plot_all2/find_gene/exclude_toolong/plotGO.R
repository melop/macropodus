setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/find_gene")
library(ggplot2)
library(clusterProfiler)

# 读取表格文件
data <- read.table("3pops.sharedHET.GO.tsv.txt", header = T, sep = "\t")

# 绘制横置条形图
ggplot(data, aes(y = reorder(Description, -pvalue), x = pvalue)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(y = "Description", x = "P-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray", size = 0.5),
        panel.grid.minor = element_line(color = "gray", size = 0.25),
        panel.grid.major.y = element_blank()) # 消除 y 轴标签和 bar 之间的网格线

