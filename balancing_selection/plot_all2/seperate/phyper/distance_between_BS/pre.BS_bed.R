setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/phyper/distance_between_BS")
library(ggplot2)
library(dplyr)
data_dir <- "/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2"  # 替换为实际的数据目录
#interval_dir <- "/data2/projects/zwang/m.hk/ROH/ROHan/runs_of_het"
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
population_order <- c("MHKmls_pop1", "MHKmls_pop2", "MHKhn", "MOPdsy", "MOPdxs_sg", "MHKjx", "MHKmlh_qns", "MOPfq_pt")
all_data$population <- factor(all_data$population, levels = population_order)
all_data <- all_data %>%
arrange(population)
# 不同种群各自betascores的top 0.1%
#horizontal_line_values <- c(MHKhk = 8.120438, MHKhn = 8.738239, MHKjx = 14.809091, MHKmlh_qns = 20.827374, MOPdsy = 11.265185, MOPfq_pt = 13.549270)
# 不同种群各自betascores的top 0.5%
horizontal_line_values <- c(MHKmls_pop1 = 30.857664, MHKmls_pop2 = 17.18587, MHKhn = 4.374229, MOPdsy = 5.915243, MOPdxs_sg = 5.836182, MHKjx = 9.462819, MHKmlh_qns = 14.154545,  MOPfq_pt = 9.55292)
all_data$hline_value <- NA
all_data$hline_value[all_data$population == "MHKmls_pop1"] <- horizontal_line_values[["MHKmls_pop1"]]
all_data$hline_value[all_data$population == "MHKmls_pop2"] <- horizontal_line_values[["MHKmls_pop2"]]
all_data$hline_value[all_data$population == "MHKhn"] <- horizontal_line_values[["MHKhn"]]
all_data$hline_value[all_data$population == "MOPdsy"] <- horizontal_line_values[["MOPdsy"]]
all_data$hline_value[all_data$population == "MOPdxs_sg"] <- horizontal_line_values[["MOPdxs_sg"]]
all_data$hline_value[all_data$population == "MHKjx"] <- horizontal_line_values[["MHKjx"]]
all_data$hline_value[all_data$population == "MHKmlh_qns"] <- horizontal_line_values[["MHKmlh_qns"]]
all_data$hline_value[all_data$population == "MOPfq_pt"] <- horizontal_line_values[["MOPfq_pt"]]
# 判断每个位点是否超过了水平线
all_data$above_threshold <- all_data$Beta1. > all_data$hline_value
all_data <- all_data %>%
mutate(window = (Position %/% 50000) * 50000)
#-------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# 定义函数用于输出每个种群的数据到文件
write_population_data <- function(population_data, population_name, output_dir) {
# 设置选项，以避免科学计数法
options(scipen = 999)

# 筛选 above_threshold 为 TRUE 的数据
filtered_data <- population_data[population_data$above_threshold, ]
filtered_data$chromosome <- gsub("chrom", "mhkscf", filtered_data$chromosome)
  
# 创建新的数据框，包含所需的三列信息
output_data <- data.frame(
  chromosome = filtered_data$chromosome,
  start = filtered_data$Position - 1,
  end = filtered_data$Position
  )    

# 设置文件名
file_name <- paste0(population_name, ".BS.bed")
file_path <- file.path(output_dir, file_name)
# 写入文件
write.table(output_data, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
# 获取输出目录
output_dir <- "/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/phyper/distance_between_BS"  # 你可以替换为其他指定的输出目录
# 为每个种群生成文件
unique_populations <- unique(all_data$population)
for (pop in unique_populations) {
pop_data <- all_data %>% filter(population == pop)
write_population_data(pop_data, pop, output_dir)
}

