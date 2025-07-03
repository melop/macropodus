#------------------------------------------------------------------------------------
# 设置工作目录
setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/long_VS_short/HIGH")

# 加载所需的包
library(ggplot2)
library(tidyr)
library(dplyr)

# 获取所有相关的文件名
files <- list.files(pattern = "genetic_load_.*\\.txt")
files <- files[grepl("FJ|PT|SZ|MHKqns|JX", files)]
datSamples <- read.table("header_samples.csv", header=F)


# 创建一个PDF文件来保存所有的图形
#pdf("allpops.refMHK.SV.Individuals_Genetic_Load_Ratio_of_ROH-HET.MOPfqpt_on_MHKjxMHKmlhqns.2.pdf", width = 10, height = 6)

# 定义种群信息
population_info <- c("MOPdsy", "MOPfq_pt", "MOPgl", "MHKhk", "MHKhn", "MHKjx", 
                     "MHKmls1", "MHKmls2", "MHKmls1", "MHKmls2", "MHKmls1", 
                     "MHKmls2", "MHKmlh_qns", "MOPdxs_sg", "MOPhn", "MOPld", 
                     "MOPdxs_sg", "MOPfq_pt", "MHKmlh_qns", "MOPyn")

times <- c(7, 1, 4, 3, 8, 4, 4, 4, 3, 1, 2, 8, 3, 2, 7, 5, 8, 2, 4, 5)
population_list <- rep(population_info, times = times)


# 读取每个文件并处理数据
#results <- list()
result_df <- data.frame(normalized_ratio = numeric(), reference_population= character(), stringsAsFactors = FALSE)
for (i in seq_along(files)) {
  file <- files[i]
#  file <- "genetic_load_PT-1.JX-4.txt"
  # 读取数据
  dat <- read.table(file, header = TRUE)
  dat <- cbind(dat, datSamples)
  
 # id_list <- unlist(str_split(gsub(".txt", "", gsub("genetic_load_", "", "genetic_load_FJ-1.JX-1.txt")), "\\."))

  parts <- strsplit(gsub("genetic_load_", "", gsub(".txt", "", file)), "\\.")[[1]]
  individual1 <- parts[1]
#  individual2 <- parts[2]
#  print(individual1)
#  print(individual2)
  # 添加种群信息到数据中
  dat$population <- population_list[1:nrow(dat)]
  
  # 计算每行的比值
  dat$individual_ratio <- dat$normalized_genetic_load_shortROH / dat$normalized_genetic_load_longROH
  
#  normalized_ratio <- dat$individual_ratio[dat$V1 == individual1] / dat$individual_ratio[dat$V1 == individual2]
  normalized_ratio <- dat$individual_ratio[dat$V1 == individual1]
    
  # 获取id_list[2]对应的population值
  reference_population <- dat$population[dat$V1 == individual1]
  #print(corresponding_population)  
  # 将当前循环计算得到的normalized_ratio和对应的population添加到结果数据框中
  result_df <- rbind(result_df, data.frame(normalized_ratio = normalized_ratio, reference_population = reference_population))
}


# 计算 p 值函数
#calculate_p_value <- function(nMean, nSD) {
#  p_value <- 1 - pnorm(1, mean = nMean, sd = nSD)
#  return(p_value)
#}
calculate_p_value <- function(nMean, nSD) {
  pnorm_value <- pnorm(1, mean = nMean, sd = nSD)
  p_value <- ifelse(pnorm_value > 0.5, 1 - pnorm_value, pnorm_value)
  return(p_value)
}

# 如果 result_df 为空，输出提示信息
if (nrow(result_df) == 0) {
  print("没有找到目标种群的数据。")
} else {
  #删除normalized_ratio为Inf的个体
  result_df <- result_df[!is.infinite(result_df$normalized_ratio), ]
  # 以参考种群进行分组，计算每个参考种群下normalized_ratio的均值、标准差和中位数等统计量
  summary_stats <- result_df %>%
    group_by(reference_population) %>%
    summarise(nMean = mean(normalized_ratio), nSD = sd(normalized_ratio, na.rm = TRUE), median_normalized_ratio = median(normalized_ratio, na.rm = TRUE), Q1 = quantile(normalized_ratio, 0.25, na.rm = TRUE)) %>%
    mutate(p_value = calculate_p_value(nMean, nSD)) %>%
    ungroup()
  
  # 将计算得到的统计量合并到result_df中
  result_df <- left_join(result_df, summary_stats, by = "reference_population")
  
  # 获取y轴的最大值
  max_y_value <- max(result_df$normalized_ratio, na.rm = TRUE)
  
  # 根据条件设置y轴的上限
  y_axis_max <- ifelse(max_y_value > 1, max_y_value + 0.05, 1.05)
  
  # 绘制箱线图，并添加p值标签
  p <- ggplot(result_df, aes(x = reference_population, y = normalized_ratio, fill = reference_population)) +
    geom_boxplot() +
    geom_text(data = summary_stats,
              aes(label = paste("p =", sprintf("%.3f", p_value)),
                  x = reference_population, y = Q1 - 0.07),
              vjust = 0, size = 8, hjust = 0.5, angle = 0, inherit.aes = FALSE) +
    labs(x = "reference_population", y = "GLoSR/GLoLR", fill = "") +
    scale_y_continuous(limits = c(NA, y_axis_max)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 24, color = "black"),
      axis.text.y = element_text(size = 24, color = "black"),
      axis.title.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24),
      plot.title = element_text(size = 24, face = "bold")
    ) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed", size = 1)
  ggsave("allpops.refMHK.NONSYN.Individuals_Genetic_Load_Ratio_of_shortROH-longROH..pdf", plot = p, width = 8, height = 6)
  print(p)
}

