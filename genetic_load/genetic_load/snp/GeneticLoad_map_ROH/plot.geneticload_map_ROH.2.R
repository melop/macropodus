# 设置工作目录
setwd("/fast3/group_crf/home/g20wangzhx36/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/HIGH/GeneticLoad_map_ROH")

# 加载所需的包
library(ggplot2)
library(lme4)

data_table <- read.csv("genetic_load_of_individual.HIGH.allpops.refMHK.csv", header=TRUE) 

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

# 定义分组
group1 <- c("MHKhn", "MOPdsy", "MOPdxs_sg", "MOPhn", "MOPyn")
group2 <- c("MHKmls_pop1", "MHKmls_pop2",  "MOPld")
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
p <- ggplot(all_individuals_data, aes(x = roh_length_Mb, y = genetic_load)) +
  geom_point(aes(color = population)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "ROH Length", y = "Genetic Load") +
  ggtitle("Linear regression analysis between genetic load and ROH length") +
  theme_minimal() +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 24)) +
  annotate("text", x = mean(range(all_individuals_data$roh_length_Mb)), y = max(all_individuals_data$genetic_load) * 0.8, label = "p=0.337", size = 8) +
  scale_color_manual(values = color_mapping)

# 修改 x 轴标签加上 "Mb"
p <- p + scale_x_continuous(labels = function(x) paste0(x, "Mb"))

print(p)
# 保存为 PDF
ggsave("Linear_regression_analysis_between_geneticload_ROHlength.HIGH.group_color.pdf", plot = p, width = 12, height = 8)

