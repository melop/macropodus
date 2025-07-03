# 设置工作目录
setwd("/fast3/group_crf/home/g20wangzhx36/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/indel/GeneticLoad_map_ROH")

# 加载所需的包
library(ggplot2)
library(lme4)

data_table <- read.csv("genetic_load_of_individual.allpops.refMHK.csv", header=TRUE) 

# 筛除 pop 为 MHKhk 和 MOPgl 的数据
data_table <- data_table[!(data_table$pop %in% c("MHKhk", "MOPgl")), ] #这两个种群没有ROH文件或者仅有一个个体存在ROH文件

#----------------------------------------------------------------------所有个体一起进行线性回归分析
# 计算个体的 ROH 长度总和函数
calculate_roh_length <- function(pop, id) {
  roh_dir <- "/fast3/group_crf/home/g20wangzhx36/m.hk/ROH/ROHan/runs_of_hom/"
  roh_file <- file.path(roh_dir, pop, paste(id, ".ROH.bed", sep = ""))
  if (file.exists(roh_file)) {
    # 检查文件大小
    if (file.info(roh_file)$size > 0) {
      roh_data <- read.table(roh_file, header = FALSE)
      if (nrow(roh_data) > 0) {
        sum(roh_data$V3 - roh_data$V2)
      } else {
        0
      }
    } else {
      0
    }
  } else {
    0
  }
}

# 收集所有有效个体的数据
all_individuals_data <- data.frame()
for (i in 1:nrow(data_table)) {
  pop <- data_table$pop[i]
  id <- data_table$id[i]
  genetic_load <- data_table$standardized_genetic_load[i]
  roh_length <- calculate_roh_length(pop, id)
  if (roh_length > 0) {
    all_individuals_data <- rbind(all_individuals_data, data.frame(genetic_load = genetic_load, roh_length = roh_length, population = pop))
  }
}

# 将 roh_length 转换为 Mb 单位
all_individuals_data$roh_length_Mb <- all_individuals_data$roh_length / 1000000

# 进行线性回归并绘图
p <- ggplot(all_individuals_data, aes(x = roh_length_Mb, y = genetic_load)) +
  geom_point(aes(color = population)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "ROH Length", y = "Genetic Load") +
  ggtitle("Linear regression analysis between genetic load and ROH length") +
  theme_minimal() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        plot.title = element_text(size = 24)) +
  annotate("text", x = mean(range(all_individuals_data$roh_length_Mb)), y = max(all_individuals_data$genetic_load) * 0.8, label = "p=0.022", size = 8)

# 修改 x 轴标签加上 "Mb"
p <- p + scale_x_continuous(labels = function(x) paste0(x, "Mb"))

print(p)
# 保存为 PDF
ggsave("Linear_regression_analysis_between_geneticload_ROHlength.pdf", plot = p, width = 12, height = 8)
#-------------------------------------------------------------------------------------混合线型模型。summary(model),将t value转换为斜率p value, p = 2*(1-pnorm(abs(t)))。0.02231283
# 计算个体的 ROH 长度总和函数
calculate_roh_length <- function(pop, id) {
  roh_dir <- "/fast3/group_crf/home/g20wangzhx36/m.hk/ROH/ROHan/runs_of_hom/"
  roh_file <- file.path(roh_dir, pop, paste(id, ".ROH.bed", sep = ""))
  if (file.exists(roh_file)) {
    # 检查文件大小
    if (file.info(roh_file)$size > 0) {
      roh_data <- read.table(roh_file, header = FALSE)
      if (nrow(roh_data) > 0) {
        sum(roh_data$V3 - roh_data$V2)
      } else {
        0
      }
    } else {
      0
    }
  } else {
    0
  }
}

# 收集所有有效个体的数据并区分不同种群颜色
all_individuals_data <- data.frame()
for (i in 1:nrow(data_table)) {
  pop <- data_table$pop[i]
  id <- data_table$id[i]
  genetic_load <- data_table$standardized_genetic_load[i]
  roh_length <- calculate_roh_length(pop, id)
  if (roh_length > 0) {
    all_individuals_data <- rbind(all_individuals_data, data.frame(genetic_load = genetic_load, roh_length = roh_length, population = pop))
  }
}

# 标准化 roh_length
all_individuals_data$roh_length_std <- scale(all_individuals_data$roh_length)

# 构建线性混合模型
model <- lmer(genetic_load ~ roh_length_std + (1 | population), data = all_individuals_data)

# 生成预测数据，确保包含 population 列
unique_populations <- unique(all_individuals_data$population)
predicted_data <- expand.grid(roh_length_std = seq(min(all_individuals_data$roh_length_std), max(all_individuals_data$roh_length_std), length.out = 100),
                              population = unique_populations)
predicted_data$genetic_load <- predict(model, newdata = predicted_data)

# 进行绘图，设置字体大小
ggplot(all_individuals_data, aes(x = roh_length_std, y = genetic_load)) +
  geom_point(aes(color = population)) +
  geom_line(data = predicted_data, aes(x = roh_length_std, y = genetic_load, color = population), linewidth = 1) +
  labs(x = "Standardized ROH Length", y = "Genetic Load") +
  ggtitle("Linear mixed model analysis between genetic load and standardized ROH length") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 20))
