setwd("/data2/projects/zwang/m.hk/ROH/ROHan")
library(ggplot2)
library(dplyr)
data=read.csv("/data2/projects/zwang/m.hk/ROH/ROHan/Froh_allDP4.csv", header = F, sep = ",")

data$Population <- factor(data$Population, levels = c("MOPdsy", "MOPld", "MOPdxs_sg", "MOPfq_pt", "MOPgl", "MOPhn", "MOPyn", "MHKhk", "MHKmlh_qns", "MHKmls", "MHKjx","MHKhn", "MSPsh", "MOCwl"))
colnames(data)[1] <- "Population" 

# 绘制箱线图
ggplot(data, aes(x = Population, y = V3, color=Population)) +
  geom_boxplot() +
  labs(x = "Population", y = "Froh") +
  theme_minimal()+
  theme(
    text = element_text(size = 30),  # 调整主要文本大小
    axis.text.x = element_text(size = 30, angle = 45, hjust = 1),  # 调整x轴文本大小和角度
    axis.text.y = element_text(size = 30),  # 调整y轴文本大小
    legend.text = element_text(size = 30),  # 调整图例文本大小
    legend.title = element_text(size = 30)  # 调整图例标题大小
  )
ggsave("Froh_allDP4.pdf", width = 16, height = 10, units = "in")
#------------------------------------------------------------------------------------------------------------------------------------------------------------
setwd("/data2/projects/zwang/m.hk/ROH/ROHan/test_MHKmls")
library(ggplot2)
library(dplyr)
data=read.csv("Froh_allDP4.MHKmls_test.csv", header = F, sep = ",")
colnames(data)[1] <- "Population" 
data <- data %>%
  filter(!grepl("MOCwl", Population))  # 将指定字符替换为您要删除的字符
data <- data %>%
  filter(!grepl("MSPsh", Population))  # 将指定字符替换为您要删除的字符
data$Population <- factor(data$Population, levels = c("MOPfq_pt", "MOPdxs_sg", "MOPld", "MOPdsy", "MOPgl", "MOPhn", "MOPyn", "MHKmls_pop2", "MHKmls_pop1", "MHKhk", "MHKmlh_qns", "MHKjx", "MHKhn"))

# 使用 aggregate 函数计算每个种群的 Froh 平均值
average_by_population <- aggregate(V3 ~ Population, data = data, FUN = mean)

print(average_by_population)

# 绘制箱线图
ggplot(data, aes(x = Population, y = V3, fill=Population)) +
  geom_boxplot() +
  labs(x = "Population", y = "Froh") +
  theme_minimal()+
  theme_classic() +
  theme(
    text = element_text(size = 45),  # 调整主要文本大小
    axis.text.x = element_text(size = 45, angle = 45, hjust = 1),  # 调整x轴文本大小和角度
    axis.text.y = element_text(size = 45),  # 调整y轴文本大小
    legend.text = element_text(size = 45),  # 调整图例文本大小
    legend.title = element_text(size = 45)  # 调整图例标题大小
  ) +
#  theme_classic() +
  coord_flip() # 翻转坐标轴
ggsave("Froh_allDP4_MHKmls_test.2.pdf", width = 16, height = 10, units = "in")
