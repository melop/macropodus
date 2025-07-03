#------------------------------------------------------------------------------------
# 设置工作目录
setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/ROH/HIGH")

# 加载所需的包
library(ggplot2)
library(tidyr)
library(dplyr)

# 获取所有相关的文件名
files <- list.files(pattern = "genetic_load_.*\\.txt")

# 创建一个PDF文件来保存所有的图形
pdf("allpops.refMHK.NONSYN.HIGH.Individuals_Genetic_Load_Ratio_of_ROH-HET.2.pdf", width = 10, height = 6)

# 定义种群信息
population_info <- c("MOPdsy", "MOPfq_pt", "MOPgl", "MHKhk", "MHKhn", "MHKjx", 
                     "MHKmls1", "MHKmls2", "MHKmls1", "MHKmls2", "MHKmls1", 
                     "MHKmls2", "MHKmlh_qns", "MOPdxs_sg", "MOPhn", "MOPld", 
                     "MOPdxs_sg", "MOPfq_pt", "MHKmlh_qns", "MOPyn")
times <- c(7, 1, 4, 3, 8, 4, 4, 4, 3, 1, 2, 8, 3, 2, 7, 5, 8, 2, 4, 5)
population_list <- rep(population_info, times = times)

# 指定参考个体的比值的索引
specified_indices <- c(8, 24, 25, 26, 27, 50, 51, 52, 75, 76, 77, 78, 79, 80) # FJ JX MHKqns PT SZ
# 存储参考个体比值的列表
reference_ratios <- rep(NA, length(specified_indices))

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
    # 从 dat 中获取对应行的比值
    if (index_in_list <= nrow(dat)) {
      reference_ratios[i] <- as.numeric(dat$individual_ratio[index_in_list])
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

# 初始化一个空的数据框，用于存储所有个体的比值
plot_data <- data.frame()

# 遍历 results 列表
for (i in seq_along(results)) {
  individual <- names(results)[i]
  for (target_pop in target_populations) {
    dat_temp <- target_data[target_data$individual == individual & target_data$population == target_pop, ]
    
    # 如果 dat_temp 不为空，则进行计算和数据处理
    if (nrow(dat_temp) > 0) {
      target_mean_roh <- mean(dat_temp$normalized_genetic_load_roh, na.rm = TRUE)
      target_mean_het <- mean(dat_temp$normalized_genetic_load_het, na.rm = TRUE)
      target_ratio <- as.numeric(target_mean_roh / (target_mean_het + 1e-10))
      
      # 创建映射表
      mapping <- list(
        "FJ" = "MOPfq_pt",
        "PT" = "MOPfq_pt",
        "JX" = "MHKjx",
        "SZ" = "MHKmlh_qns",
        "MHKqns" = "MHKmlh_qns"
      )
      
      # 提取种群名并应用映射
      population_mapped <- sapply(individual, function(individual) {
        base_name <- strsplit(individual, "-")[[1]][1]
        if (base_name %in% names(mapping)) {
          return(mapping[[base_name]])
        } else {
          return(base_name)
        }
      })
      
      # 获取当前索引 i 对应的指定索引列表中的值
      reference_ratio <- as.numeric(reference_ratios[i])
      
      # 如果 reference_ratio 存在且为有效数字，则计算目标种群比值除以个体比值
      if (!is.na(reference_ratio) && reference_ratio != 0) {
#        normalized_ratio <- target_ratio / reference_ratio
        normalized_ratio <- reference_ratio / target_ratio
        
        # 将结果添加到数据框
        plot_data <- rbind(plot_data, data.frame(Population = population_mapped, Individual = individual, Target_pop = target_pop, Ratio = normalized_ratio))
      }
    }
  }
}

# 计算 p 值函数
calculate_p_value <- function(nMean, nSD) {
  p_value <- 1 - pnorm(1, mean = nMean, sd = nSD)
  return(p_value)
}
# 如果 plot_data 为空，输出提示信息
if (nrow(plot_data) == 0) {
  print("没有找到目标种群的数据。")
} else {
  # 计算每个种群和目标种群组合的 p 值和中位数
  summary_stats <- plot_data %>%
    group_by(Population, Target_pop) %>%
    summarise(nMean = mean(Ratio), nSD = sd(Ratio, na.rm = TRUE), median_ratio = median(Ratio, na.rm = TRUE), Q1 = quantile(Ratio, 0.25, na.rm = TRUE)) %>%
    mutate(p_value = calculate_p_value(nMean, nSD)) %>%
    ungroup()
  
  # 合并均值、标准差和 p 值到 plot_data 中
  plot_data <- left_join(plot_data, summary_stats, by = c("Population", "Target_pop"))
  
  # 获取 y 轴的最大值
  max_y_value <- max(plot_data$Ratio, na.rm = TRUE)
  
  # 根据条件设置 y 轴的上限
  y_axis_max <- ifelse(max_y_value > 1, max_y_value + 0.05, 1.05)
  
  # 绘制箱线图，并添加 p 值标签
  p <- ggplot(plot_data, aes(x = Population, y = Ratio, fill = Target_pop)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +  # 调整箱线图的位置，使其水平并列
    geom_text(data = summary_stats,
              aes(label = paste("p =", sprintf("%.3f", p_value)),
                  x = Population, y = Q1-0.08,
                  group = interaction(Population, Target_pop)),
              position = position_dodge(width = 0.8),  # 根据箱线图的位置进行调整
              vjust = 0, size = 6, hjust = 0.5, angle = 0, inherit.aes = FALSE) +
    labs(x = "Target Population", y = "(tGLroh/tGLhet)/(rmGLroh/rmGLhet)", fill = "Reference Population") +
    #上面的代码把"MOPyn", "MOPhn", "MHKhn", "MHKhk"当作reference populations，把MHKmlh_qns,MHKjx和MOPfq_pt当作target populations，
    #实际上应该反过来，因此仅在最后绘图这一行代码将二者的名称对调即可。
    # 使用 theme_minimal()，并在其中调整字体大小
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 16),  # x轴刻度标签的字体大小
      axis.text.y = element_text(size = 16),  # y轴刻度标签的字体大小
      axis.title.x = element_text(size = 20),  # x轴标题的字体大小
      axis.title.y = element_text(size = 20),  # y轴标题的字体大小
      legend.text = element_text(size = 16),  # 图例文本的字体大小
      legend.title = element_text(size = 16),  # 图例标题的字体大小
      plot.title = element_text(size = 24, face = "bold")  # 图形标题的字体大小和加粗样式
    ) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed", size = 1)  # 添加红色水平线
  
  print(p)  # 直接打印绘制的图形
}

# 关闭 PDF 文件
dev.off()

