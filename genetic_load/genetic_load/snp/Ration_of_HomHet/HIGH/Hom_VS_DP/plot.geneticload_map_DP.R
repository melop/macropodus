#------------------------------------------------------------------
# 设置工作目录
setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/Ration_of_HomHet/HIGH/Hom_VS_DP")

# 加载所需的包
library(ggplot2)
library(lme4)

# 读取表格数据
table1 <- read.table("genotype_proportions.HIGH.txt", header = T)
samples_ID <- read.csv("header_samples.csv", header = F, col.names = "id")
table1_ID <- cbind(table1, samples_ID)
table2 <- read.csv("all_MHK_MOP.depth.csv", header = F, col.names = c("rid","DP","id","pop"))


# 合并两个表格，基于个体名进行匹配
data_table <- merge(table1_ID, table2, by = "id", all.x = FALSE)


#----------------------------------------------------------------------所有个体一起进行线性回归分析

# 收集所有有效个体的数据
all_individuals_data <- data.frame()
for (i in 1:nrow(data_table)) {
  pop <- data_table$pop[i]
  id <- data_table$id[i]
  DP <- data_table$DP[i]
  HomPP <- data_table$Proportion[i]
  all_individuals_data <- rbind(all_individuals_data, data.frame(HomPP = HomPP, id = id, population = pop, DP = DP))
}


# 定义分组
group1 <- c("MHKhn", "MHKhk", "MOPdsy", "MOPdxs_sg", "MOPgl","MOPhn", "MOPyn")
group2 <- c("MHKgbc", "MHKlw",  "MOPld")
group3 <- c("MHKjx", "MHKmlh_qns", "MOPfq_pt")

# 为每个组选择一组相近的颜色，这里使用不同的色系来区分组
color_group1 <- colorRampPalette(c("blue", "cyan"))(length(group1))
color_group3 <- colorRampPalette(c("red", "orange"))(length(group2))
color_group2 <- colorRampPalette(c("green", "yellowgreen"))(length(group3))

# 创建颜色映射
color_mapping <- c(setNames(color_group1, group1),
                   setNames(color_group2, group2),
                   setNames(color_group3, group3))

# 将 population 列转换为有序因子，按照分组顺序设置水平
all_individuals_data$population <- factor(all_individuals_data$population, 
                                          levels = c(group1, group2, group3))

# 进行线性回归并绘图
p <- ggplot(all_individuals_data, aes(x = DP, y = HomPP)) +
  geom_point(aes(color = population)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "DP", y = "Homozygosity Propotion") +
  ggtitle("Linear regression analysis between Homozygosity Propotion and DP") +
  theme_minimal() +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 24)) +
  annotate("text", x = mean(range(all_individuals_data$DP)), y = max(all_individuals_data$HomPP) * 0.8, label = "p=0.1076", size = 8) +
  scale_color_manual(values = color_mapping)


print(p)
# 保存为 PDF
ggsave("Linear_regression_analysis_between_HomPP_DP.group_color.snp.HIGH.pdf", plot = p, width = 12, height = 8)

#-------------------------------------------------------------------------------------混合线型模型。summary(model),将t value转换为斜率p value, p = 2*(1-pnorm(abs(t)))。0.1076

# 收集所有有效个体的数据并区分不同种群颜色
all_individuals_data <- data.frame()
for (i in 1:nrow(data_table)) {
  pop <- data_table$pop[i]
  id <- data_table$id[i]
  DP <- data_table$DP[i]
  HomPP <- data_table$Proportion[i]
  all_individuals_data <- rbind(all_individuals_data, data.frame(HomPP = HomPP, id = id, population = pop, DP = DP))
}

# 标准化 roh_length
all_individuals_data$DP_std <- scale(all_individuals_data$DP)

# 构建线性混合模型
model <- lmer(HomPP ~ DP_std + (1 | population), data = all_individuals_data)

# 生成预测数据，确保包含 population 列
unique_populations <- unique(all_individuals_data$population)
predicted_data <- expand.grid(DP_std = seq(min(all_individuals_data$DP_std), max(all_individuals_data$DP_std), length.out = 100),
                              population = unique_populations)
predicted_data$HomPP <- predict(model, newdata = predicted_data)

# 进行绘图，设置字体大小
ggplot(all_individuals_data, aes(x = DP_std, y = HomPP)) +
  geom_point(aes(color = population)) +
  geom_line(data = predicted_data, aes(x = DP_std, y = HomPP, color = population), size = 1) +
  labs(x = "Standardized DP", y = "Homozygosity Propotion") +
  ggtitle("Linear mixed model analysis between Homozygosity Propotion and DP") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 20))

