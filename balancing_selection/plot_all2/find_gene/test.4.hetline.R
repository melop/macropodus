setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/find_gene")

library(ggplot2)
library(dplyr)

data_dir <- "/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2"  # 替换为实际的数据目录
interval_dir <- "/data2/projects/zwang/m.hk/ROH/ROHan/runs_of_het"

# 读取文件数据
read_chromosome_data <- function(file_path) {
  data <- read.table(file_path, header = TRUE, na.strings = "NA")
  # 提取种群和染色体编号
  file_name <- basename(file_path)
  pop_chr <- gsub("\\.betascores\\.txt", "", file_name)
  pop_chr <- strsplit(pop_chr, "\\.")[[1]]
  pop_id <- as.character(pop_chr[1])
  chr_id <- as.character(pop_chr[2])
  data$chromosome <- chr_id
  data$population <- pop_id
  return(data)
}

# 读取所有文件数据
file_paths <- list.files(data_dir, pattern = "\\.betascores\\.txt$", full.names = TRUE)
all_data <- do.call(rbind, lapply(file_paths, read_chromosome_data))
# 筛选出三个种群的数据
part_data <- all_data[all_data$population %in% c("MHKjx", "MHKmlh_qns", "MOPfq_pt"), ]

#population_order <- c("MHKjx", "MHKmlh_qns", "MOPfq_pt")
#all_data$population <- factor(all_data$population, levels = population_order)
#all_data <- all_data %>%
#  arrange(population)

# 根据种群名称读取区间数据，并为不同个体进行标记
read_interval_data_for_population <- function(pop_id) {
  # 获取对应的区间文件
  interval_files <- list.files(file.path(interval_dir, pop_id), pattern = "\\.het\\.bed$", full.names = TRUE)  
  interval_data_list <- lapply(interval_files, function(file) {
    # 读取区间数据
    interval_data <- read.table(file, header = FALSE, col.names = c("chromosome", "start", "end", "het"))
    
    # 提取个体ID，假设个体ID是文件名的一部分
    individual_id <- tools::file_path_sans_ext(basename(file))  # 获取文件名（不带扩展名）
    
    # 将种群和个体信息添加到数据中
    interval_data$population <- pop_id
    interval_data$individual <- individual_id  # 添加个体标识
    
    return(interval_data)
  })
  
  # 合并所有个体的数据
  return(do.call(rbind, interval_data_list))
}

# 为所有种群读取区间文件并合并
unique_populations <- unique(part_data$population)
interval_data_all <- do.call(rbind, lapply(unique_populations, read_interval_data_for_population))
interval_data_all$chromosome <- gsub("mhkscf", "chrom", interval_data_all$chromosome)

# 不同种群各自betascores的top 0.1%
#horizontal_line_values <- c(MHKhk = 8.120438, MHKhn = 8.738239, MHKjx = 14.809091, MHKmlh_qns = 20.827374, MOPdsy = 11.265185, MOPfq_pt = 13.549270)

# 不同种群各自betascores的top 0.5%
horizontal_line_values <- c(MHKjx = 9.462819, MHKmlh_qns = 14.154545,  MOPfq_pt = 9.55292)

part_data$hline_value <- NA
part_data$hline_value[part_data$population == "MHKjx"] <- horizontal_line_values[["MHKjx"]]
part_data$hline_value[part_data$population == "MHKmlh_qns"] <- horizontal_line_values[["MHKmlh_qns"]]
part_data$hline_value[part_data$population == "MOPfq_pt"] <- horizontal_line_values[["MOPfq_pt"]]

#------------------------------------------------------------------------
# 计算每个种群的个体数量
individual_counts <- interval_data_all %>%
  group_by(population) %>%
  summarise(num_individuals = n_distinct(individual))  # 使用个体标识符作为标识

# 为区间横线设置偏移量，使得每个个体的横线位于不同的y水平线
interval_data_all <- interval_data_all %>%
  group_by(population) %>%
  mutate(
    individual_rank = dense_rank(individual),  # 使用密集排名
    base_y_position = horizontal_line_values[population] - 2,
    y_position = base_y_position - (individual_rank - 1) * 0.5
  ) %>%
  ungroup()

# 确保y_position不会低于阈值线减去一定偏移量
interval_data_all$y_position <- pmax(interval_data_all$y_position, horizontal_line_values[interval_data_all$population] - 5)

#------------------------------------------------------------------------

# 判断每个位点是否超过了水平线
part_data$above_threshold <- part_data$Beta1. > part_data$hline_value

# 过滤出大于阈值的数据
filtered_data <- part_data %>% filter(above_threshold)

# 添加窗口分组
filtered_data <- filtered_data %>%
  mutate(window = (Position %/% 50000) * 50000)

# 统计每个窗口内超过阈值的种群数量
threshold_summary <- filtered_data %>%
  group_by(chromosome, window) %>%
  summarise(count_above = sum(above_threshold), 
            num_populations_above = n_distinct(population[above_threshold])) %>%
  ungroup()

# 合并统计结果
filtered_data <- left_join(filtered_data, threshold_summary, by = c("chromosome", "window"))

# 标记highlight
filtered_data$highlight <- filtered_data$num_populations_above >= 2

# 筛选出MHKjx, MHKmlh_qns, MOPfq_pt三个种群的第6号染色体
filtered_data_subset <- filtered_data %>%
  filter(population %in% c("MHKjx", "MHKmlh_qns", "MOPfq_pt") & chromosome == "chrom_6")

# 筛选区间数据
interval_data_subset <- interval_data_all %>%
  filter(population %in% c("MHKjx", "MHKmlh_qns", "MOPfq_pt") & chromosome == "chrom_6")

# 绘制曼哈顿图（只包含三个种群的第6号染色体）
p_subset <- ggplot(filtered_data_subset, aes(x = Position, y = Beta1.)) +
  geom_point(aes(color = factor(chromosome)), alpha = 0.6, size = 0.3) +
  geom_point(data = filtered_data_subset[filtered_data_subset$highlight == TRUE, ], aes(x = Position, y = Beta1.), color = "black", size = 1) +
  facet_grid(population ~ chromosome, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = "Position",
    y = "Beta Value",
    title = "Manhattan Plot for Chromosome 6 of MHKjx, MHKmlh_qns, and MOPfq_pt"
  )

# 添加水平线
p_subset <- p_subset + geom_hline(data = filtered_data_subset[filtered_data_subset$population == "MHKjx", ], aes(yintercept = horizontal_line_values[["MHKjx"]]), linetype = "dashed", color = "red") +
  geom_hline(data = filtered_data_subset[filtered_data_subset$population == "MHKmlh_qns", ], aes(yintercept = horizontal_line_values[["MHKmlh_qns"]]), linetype = "dashed", color = "red") +
  geom_hline(data = filtered_data_subset[filtered_data_subset$population == "MOPfq_pt", ], aes(yintercept = horizontal_line_values[["MOPfq_pt"]]), linetype = "dashed", color = "red")

# 添加区间横线
p_subset <- p_subset + geom_segment(data = interval_data_subset, aes(
  x = start,
  xend = end,
  y = y_position,
  yend = y_position
), color = "blue", size = 1.5)

# 保存为PDF
ggsave("manhattan_plot_subset_chrom6.2.pdf", plot = p_subset, width = 10, height = 6)
print(p_subset)


