# 设置工作目录
setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/indel/ROH")

# 加载所需的包
library(ggplot2)
library(tidyr)

# 获取所有相关的文件名
files <- list.files(pattern = "genetic_load_.*\\.txt")

# 创建一个PDF文件来保存所有的图形
pdf("allpops.refMHK.SV.Individuals_Genetic_Load_Ratio_of_ROH-HET.pdf", width = 10, height = 6)

# 定义种群信息
population_info <- c("MOPdsy", "MOPfq_pt", "MOPgl", "MHKhk", "MHKhn", "MHKjx", 
                     "MHKmls1", "MHKmls2", "MHKmls1", "MHKmls2", "MHKmls1", 
                     "MHKmls2", "MHKmlh_qns", "MOPdxs_sg", "MOPhn", "MOPld", 
                     "MOPdxs_sg", "MOPfq_pt", "MHKmlh_qns", "MOPyn")
times <- c(7, 1, 4, 3, 8, 4, 4, 4, 3, 1, 2, 8, 3, 2, 7, 5, 8, 2, 4, 5)
population_list <- rep(population_info, times = times)

# 指定参考个体的比值的索引
specified_indices <- c(8, 24, 25, 26, 27, 50, 51, 52, 75, 76, 77, 78, 79, 80) #FJ JX MHKqns PT SZ
# 存储参考个体比值的列表
reference_ratios <- list()

# 读取每个文件并处理数据
results <- list()

for (i in seq_along(files)) {
  file <- files[i]
  # 读取数据
  dat <- read.table(file, header = TRUE)
  
  # 提取种群名和个体编号
  filename_parts <- strsplit(basename(file), "[-_.]")[[1]]
  population <- filename_parts[3]  # 提取种群名，比如 "JX"
  individual <- paste(filename_parts[3], filename_parts[4], sep = "-")  # 提取完整个体名，比如 "JX-1"
  
  # 添加种群信息到数据中
  dat$population <- population_list[1:nrow(dat)]
  
  # 计算每行的比值
  dat$individual_ratio <- dat$normalized_genetic_load_roh / dat$normalized_genetic_load_het
  
  # 获取当前索引 i 对应的指定索引列表中的值
  if (i <= length(specified_indices)) {
    index_in_list <- specified_indices[i]
    # 从 dat_temp 中获取对应行的比值
    if (index_in_list <= nrow(dat)) {
      reference_ratios[[i]] <- dat$individual_ratio[index_in_list]
    }
  }
  # 添加个体信息
  dat$individual <- individual
  
  # 保存结果到列表中
  results[[individual]] <- dat
}

# 合并目标种群的数据
# 目标种群列表
target_populations <- c("MOPyn", "MOPhn", "MHKhn", "MHKhk")

# 使用 lapply 函数遍历每个数据框，并筛选出目标种群的行
target_data <- do.call(rbind, lapply(results, function(dat) {
  # 筛选出属于目标种群的行
  target_rows <- dat[dat$population %in% target_populations, ]
  
  # 如果目标种群行存在，返回这些行
  if (nrow(target_rows) > 0) {
    return(target_rows)
  }
}))

# 如果 target_data 为空，输出提示信息
if (nrow(target_data) == 0) {
  print("没有找到目标种群的数据。")
}


# 初始化一个空的数据框，用于存储所有个体的比值
plot_data <- data.frame()

# 遍历 results 列表
for (i in seq_along(results)) {
  individual <- names(results)[i]
  # 获取个体数据
  dat_temp <- target_data[target_data$individual == individual, ]
  
  # 计算目标种群的genetic_load_ROH和genetic_load_HET平均值的比值
  target_mean_roh <- mean(dat_temp$normalized_genetic_load_roh)
  target_mean_het <- mean(dat_temp$normalized_genetic_load_het)
  target_ratio <- target_mean_roh / (target_mean_het + 1e-10)
  
  # 创建映射表
  mapping <- list(
    "FJ" = "MOPfq_pt",
    "PT" = "MOPfq_pt",
    "JX" = "MHKjx",
    "SZ" = "MHKmlh_qns",
    "MHKqns" = "MHKmlh_qns"
  )
  # 提取种群名并应用映射
  population <- sapply(individual, function(individual) {
    # 提取种群名
    base_name <- strsplit(individual, "-")[[1]][1]
    
    # 根据映射表转换种群名
    if (base_name %in% names(mapping)) {
      return(mapping[[base_name]])
    } else {
      return(base_name)  # 如果不在映射表中，则返回原种群名
    }
  })
  
  # 获取当前索引 i 对应的指定索引列表中的值
  reference_ratio <- as.numeric(reference_ratios[i])
  
  # 计算目标种群比值除以个体比值
  normalized_ratio <- target_ratio / reference_ratio
  
  # 将结果添加到数据框
  plot_data <- rbind(plot_data, data.frame(Population = population, Individual = individual, Ratio = normalized_ratio))
}

# 绘制按种群分类的箱线图
p <- ggplot(plot_data, aes(x = Population, y = Ratio, fill = Population)) +
  geom_boxplot() +
  labs(title = paste("Normalized Genetic Load Ratio by Population (Target Populations:", paste(target_populations, collapse = ", "), ")"),
       x = "Population",
       y = "Normalized Genetic Load Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# 显示图形
print(p)

# 保存图形为PDF
ggsave("genetic_load_ratio_boxplot.pdf", plot = p, width = 10, height = 6)