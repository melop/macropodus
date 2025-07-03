setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2")

library(ggplot2)
library(dplyr)

data_dir <- "/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2"  # 替换为实际的数据目录

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

# 不同种群各自betascores的top 0.1%
#horizontal_line_values <- c(MHKhk = 8.120438, MHKhn = 8.738239, MHKjx = 14.809091, MHKmlh_qns = 20.827374, MOPdsy = 11.265185, MOPfq_pt = 13.549270)

# 不同种群各自betascores的top 0.5%
horizontal_line_values <- c(MHKmls_pop1 = 30.857664, MHKmls_pop2 = 17.18587, MHKhn = 4.374229, MHKjx = 9.462819, MHKmlh_qns = 14.154545, MOPdsy = 5.915243, MOPfq_pt = 9.55292, MOPdxs_sg = 5.836182)

all_data$hline_value <- NA
all_data$hline_value[all_data$population == "MHKmls_pop1"] <- horizontal_line_values[["MHKmls_pop1"]]
all_data$hline_value[all_data$population == "MHKmls_pop2"] <- horizontal_line_values[["MHKmls_pop2"]]
all_data$hline_value[all_data$population == "MHKhn"] <- horizontal_line_values[["MHKhn"]]
all_data$hline_value[all_data$population == "MHKjx"] <- horizontal_line_values[["MHKjx"]]
all_data$hline_value[all_data$population == "MHKmlh_qns"] <- horizontal_line_values[["MHKmlh_qns"]]
all_data$hline_value[all_data$population == "MOPdsy"] <- horizontal_line_values[["MOPdsy"]]
all_data$hline_value[all_data$population == "MOPfq_pt"] <- horizontal_line_values[["MOPfq_pt"]]
all_data$hline_value[all_data$population == "MOPdxs_sg"] <- horizontal_line_values[["MOPdxs_sg"]]

# 判断每个位点是否超过了水平线
all_data$above_threshold <- all_data$Beta1. > all_data$hline_value

# 按50000 bp窗口进行分组
all_data <- all_data %>%
  mutate(window = (Position %/% 50000) * 50000)


# 统计每个窗口内超过阈值的种群数量
threshold_summary <- all_data %>%
  group_by(chromosome, window) %>%
  summarise(count_above = sum(above_threshold), 
            num_populations_above = n_distinct(population[above_threshold])) %>%
  ungroup()


# 将统计结果合并回原数据框
all_data <- left_join(all_data, threshold_summary, by = c("chromosome", "window"))
#-----------------------------------------------------------------------------------------上面处理好表格，下面可以进行各种统计
# 如果 count_above >= 2 且 num_populations_above >= 2，则标记为 TRUE，否则为 FALSE
all_data$highlight <- all_data$num_populations_above >= 2 & all_data$above_threshold

# 指定种群的排列顺序
population_order <- c("MHKmls_pop1", "MHKmls_pop2","MHKhn", "MOPdsy", "MOPdxs_sg", "MHKjx", "MHKmlh_qns", "MOPfq_pt")  # 替换为实际的种群ID按所需顺序排列
all_data$population <- factor(all_data$population, levels = population_order)

# 绘制曼哈顿图，标记符合条件的位点为黑色
p <- ggplot(all_data, aes(x = Position, y = Beta1.)) +
  geom_point(aes(color = factor(chromosome)), alpha = 0.6, size = 0.5) +  # 调整点的大小和透明度
  geom_point(data = all_data[all_data$highlight == TRUE, ], aes(x = Position, y = Beta1.), color = "black", size = 1) +  # 标记满足条件的位点
  facet_grid(population ~ chromosome, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = "Position",
    y = "Beta Value",
    title = "Manhattan Plots of Beta Values Across Chromosomes for Each Population"
  )

# 添加每个种群的水平线
p <- p + geom_hline(data = all_data[all_data$population == "MHKmls_pop1", ], aes(yintercept = horizontal_line_values[["MHKmls_pop1"]]), linetype = "dashed", color = "red") +
  geom_hline(data = all_data[all_data$population == "MHKmls_pop2", ], aes(yintercept = horizontal_line_values[["MHKmls_pop2"]]), linetype = "dashed", color = "red") +
  geom_hline(data = all_data[all_data$population == "MHKhn", ], aes(yintercept = horizontal_line_values[["MHKhn"]]), linetype = "dashed", color = "red") +
  geom_hline(data = all_data[all_data$population == "MHKjx", ], aes(yintercept = horizontal_line_values[["MHKjx"]]), linetype = "dashed", color = "red") +
  geom_hline(data = all_data[all_data$population == "MHKmlh_qns", ], aes(yintercept = horizontal_line_values[["MHKmlh_qns"]]), linetype = "dashed", color = "red") +
  geom_hline(data = all_data[all_data$population == "MOPdsy", ], aes(yintercept = horizontal_line_values[["MOPdsy"]]), linetype = "dashed", color = "red") +
  geom_hline(data = all_data[all_data$population == "MOPfq_pt", ], aes(yintercept = horizontal_line_values[["MOPfq_pt"]]), linetype = "dashed", color = "red") +
  geom_hline(data = all_data[all_data$population == "MOPdxs_sg", ], aes(yintercept = horizontal_line_values[["MOPdxs_sg"]]), linetype = "dashed", color = "red")
# 保存曼哈顿图为PDF文件
ggsave("manhattan_plot.pdf", plot = p, width = 10, height = 6)
print(p)
