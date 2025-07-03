setwd("/data2/projects/zwang/m.hk/ROH/ROHan/50kb_ROH")
# 读取数据
froh_data <- read.csv("Froh_allDP4.MHKmls_test.csv", header = F, col.names = c("population", "individual", "froh"))


# 按种群分组计算平均值
population_stats <- aggregate(
  froh ~ population,  # 替换为你的列名（格式："值列 ~ 分组列"）
  data = froh_data,
  FUN = mean,
  na.rm = TRUE
)

# 重命名列（可选）
colnames(population_stats) <- c("population", "mean_froh")

# 保存结果
write.csv(population_stats, "allpops.Froh_50kb.csv", row.names = FALSE)
