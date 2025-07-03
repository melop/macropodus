setwd("/fast3/group_crf/home/g20wangzhx36/m.hk/Structure/Pixy/output/dxy_bed/dxy_map_ROH")

library(tidyverse)
# 查找当前目录下所有子目录中的后缀为meandxy.HET.HOM.txt的文件
file_list <- list.files(path = ".", pattern = "meandxy.HET.HOM.txt", recursive = TRUE, full.names = TRUE)

# 读取所有符合条件的文件并合并成一个数据框，同时添加种群组合列
data <- map_df(file_list, function(x) {
  file_data <- read.table(x, header = FALSE, sep = "\t")
  colnames(file_data) <- c("Individual", "Region_Type", "meandxy")
  population_combination <- sub("(.*)\\.meandxy\\.HET\\.HOM\\.txt", "\\1", basename(x))
  file_data <- file_data %>%
    mutate(Population_Combination = population_combination)
  return(file_data)
})
# 数据预处理
# 添加一个新列用于区分HET和HOM区域（提取Region_Type中的关键信息）
data <- data %>%
  mutate(Region = ifelse(grepl("HET", Region_Type), "HET", "HOM"))

# 绘制箱线图
ggplot(data, aes(x = Region, y = meandxy, fill = Region)) +
  geom_boxplot() +
  facet_wrap(~Population_Combination) +
  labs(x = "Region", y = "Mean dxy", title = "Mean dxy by Region for Each Population Combination") +
  theme_minimal()

# 进行显著性检验并标记p值
p_values <- data %>%
  group_by(Population_Combination) %>%
  summarise(p_value = wilcox.test(meandxy ~ Region)$p.value)

