setwd("/fast3/group_crf/home/g20wangzhx36/m.hk/Structure/Pixy/output/fst_bed")

# 加载所需的包
library(ggplot2)
library(dplyr)
library(colorspace)

# 定义文件列表，这里可根据实际情况添加或修改文件名称
file_list <- c("MOPdsy.MHKhn", "MOPdsy.MHKmlh_qns", "MOPfq_pt.MHKhn", "MOPfq_pt.MHKmlh_qns")

# 设置文件后缀，假设数据文件的实际完整名称是类似 "MOPdsy.MHKhn.dxy.bed" 这种格式，根据实际情况修改
file_suffix <- ".fst.bed"

# 创建一个pdf文件用于保存所有的图
pdf("presentive_MHK_MOP.fst.manhattan.pdf", width = 18, height = 6)

for (file_base_name in file_list) {
  # 构建完整的文件名
  file_name <- paste0(file_base_name, file_suffix)
  
  # 读取数据文件
  data <- read.table(file_name, header = FALSE, sep = "\t", 
                     col.names = c("chromosome", "start", "end", "fst"))
  data$chromosome <- as.numeric(gsub("mhkscf_", "", data$chromosome))
  
  # 筛选出染色体编号小于24的窗口数据
  filtered_data <- data[data$chromosome < 24, ]
  
  # 进一步筛选，去除fst为NA的行
  filtered_data <- filtered_data[!is.na(filtered_data$fst), ]
  
  # 计算窗口中点坐标
  filtered_data <- filtered_data %>%
    mutate(midpoint = (start + end) / 2)
  
  # 对染色体列进行因子化处理
  chromosome_levels <- unique(filtered_data$chromosome)
  filtered_data$chromosome <- factor(filtered_data$chromosome, levels = chromosome_levels)
  
  # 计算每个染色体区域的范围
  chromosome_ranges <- filtered_data %>%
    group_by(chromosome) %>%
    summarize(min_start = min(start), max_end = max(end)) %>%
    mutate(chromosome_start = cumsum(c(0, lag(max_end, default = 0)[-1])),
           chromosome_end = chromosome_start + (max_end - min_start))
  
  # 将窗口中点坐标映射到对应染色体区域的相对坐标位置
  filtered_data <- filtered_data %>%
    left_join(chromosome_ranges, by = "chromosome") %>%
    mutate(adjusted_midpoint = chromosome_start + (midpoint - min_start))
  
  # 找到当前文件fst的top1%的值
  top_1_percent_value <- quantile(filtered_data$fst, 0.99, na.rm = TRUE)
  
  # 根据fst与top1%值的比较，设置颜色列
  filtered_data <- filtered_data %>%
    mutate(point_color = ifelse(fst >= top_1_percent_value, as.character(chromosome), "gray"))
  # 将小于阈值的点颜色设置为灰色
  #filtered_data$point_color[filtered_data$point_color == "gray"] <- "gray"
  # 获取不同染色体编号对应的rainbow颜色向量（这里假设chromosome列是整数编号）
  chromosome_levels <- unique(filtered_data$chromosome)
  rainbow_colors <- rainbow(length(chromosome_levels), alpha = 0.8, start = 0, end = 1)
  
  # 构建颜色映射的名称和颜色对应关系的数据框（用于scale_color_manual）
  color_mapping <- data.frame(chromosome = as.character(chromosome_levels), color = rainbow_colors)
  
  # 绘制曼哈顿图
  p <- ggplot(filtered_data, aes(x = adjusted_midpoint, y = fst, color = point_color)) +
    geom_point(size = 1, alpha = 0.8) +
    labs(x = "Chromosome", y = "fst", title = file_base_name) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "none",
          text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x = element_blank()) +  
    scale_x_continuous(breaks = chromosome_ranges$chromosome_start + (chromosome_ranges$max_end - chromosome_ranges$min_start) / 2,
                       labels = chromosome_ranges$chromosome) +
    geom_hline(yintercept = top_1_percent_value, color = "red", linetype = "dashed") +
    scale_color_manual(values = setNames(rainbow_colors, as.character(chromosome_levels)))
  
  print(p)
}

# 关闭pdf文件设备，结束保存操作
dev.off()


