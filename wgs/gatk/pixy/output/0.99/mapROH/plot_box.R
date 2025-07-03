setwd("/fast3/group_crf/home/g20wangzhx36/m.hk/Structure/Pixy/output/0.99/mapROH")

library(tidyverse)

# 查找当前目录下所有子目录中后缀为het.dxy.bed、hom.dxy.bed、het.bed或hom.bed的文件
file_list <- list.files(path = ".", pattern = "het.dxy.bed|hom.dxy.bed|het.bed|hom.bed", recursive = TRUE, full.names = TRUE)

# 创建一个空的数据框来存储个体名、种群名和各文件行数
result_df <- data.frame(Individual_Name = character(),
                        Population_Name = character(),
                        Het_dxy_bed_Rows = integer(),
                        Hom_dxy_bed_Rows = integer(),
                        Het_bed_Rows = integer(),
                        Hom_bed_Rows = integer(),
                        stringsAsFactors = FALSE)

# 遍历文件列表
for (i in 1:length(file_list)) {
  file_path <- file_list[i]
  # 获取文件名中的个体名，考虑四种后缀情况
  individual_name <- sub("(.*)\\.(het|hom)\\.(dxy|)\\.bed", "\\1", basename(file_path))
  # 获取文件所在的子目录名作为种群名
  population_name <- basename(dirname(file_path))
  
  if (grepl("het.dxy.bed", file_path)) {
    het_dxy_bed_rows <- nrow(read.table(file_path, header = FALSE))
    # 查找对应的hom.dxy.bed文件
    hom_dxy_bed_file_path <- file_list[grepl(paste0(individual_name, "\\.hom.dxy.bed"), file_list)]
    if (length(hom_dxy_bed_file_path) == 1) {
      hom_dxy_bed_rows <- nrow(read.table(hom_dxy_bed_file_path, header = FALSE))
      # 查找对应的het.bed文件
      het_bed_file_path <- file_list[grepl(paste0(individual_name, "\\.het.bed"), file_list)]
      if (length(het_bed_file_path) == 1) {
        het_bed_rows <- nrow(read.table(het_bed_file_path, header = FALSE))
        # 查找对应的hom.bed文件
        hom_bed_file_path <- file_list[grepl(paste0(individual_name, "\\.hom.bed"), file_list)]
        if (length(hom_bed_file_path) == 1) {
          hom_bed_rows <- nrow(read.table(hom_bed_file_path, header = FALSE))
          
          # 将个体名、种群名和各文件行数添加到结果数据框中
          new_row <- data.frame(Individual_Name = individual_name,
                                Population_Name = population_name,
                                Het_dxy_bed_Rows = het_dxy_bed_rows,
                                Hom_dxy_bed_Rows = hom_dxy_bed_rows,
                                Het_bed_Rows = het_bed_rows,
                                Hom_bed_Rows = hom_bed_rows,
                                stringsAsFactors = FALSE)
          result_df <- rbind(result_df, new_row)
        }
      }
    }
  }
}
# 查看结果
result_df

library(dplyr)
# 按种群合并个体数据
merge_data_by_population <- function(df) {
  merged_df <- df %>%
    group_by(Population_Name) %>%
    summarise(Het_dxy_bed_Rows = sum(Het_dxy_bed_Rows),
              Hom_dxy_bed_Rows = sum(Hom_dxy_bed_Rows),
              Het_bed_Rows = sum(Het_bed_Rows),
              Hom_bed_Rows = sum(Hom_bed_Rows))
  return(merged_df)
}

# 构建列联表的函数
build_contingency_table <- function(population_data) {
  contingency_table <- matrix(c(population_data$Het_dxy_bed_Rows,
                                population_data$Hom_dxy_bed_Rows,
                                population_data$Het_bed_Rows,
                                population_data$Hom_bed_Rows),
                              nrow = 2, byrow = TRUE)
  colnames(contingency_table) <- c("dxy", "total")
  rownames(contingency_table) <- c("Het", "Hom")
  return(contingency_table)
}
# 对种群数据进行卡方检验的函数
chi_square_test_per_population <- function(population_data) {
  contingency_table <- build_contingency_table(population_data)
  chi_square_result <- chisq.test(contingency_table)
  return(chi_square_result$p.value)
}

# 先按种群合并个体数据
merged_result_df <- merge_data_by_population(result_df)

# 按照种群分组并进行卡方检验
p_values <- merged_result_df %>%
  group_by(Population_Name) %>%
  do(p_value = chi_square_test_per_population(.))

# 将p-value添加到merged_result_df对应的行
merged_result_df <- merged_result_df %>%
  left_join(p_values, by = "Population_Name")

# 重命名列名以便更清晰
colnames(merged_result_df)[colnames(merged_result_df) == "p_value"] <- "Chi_Square_P_Value"


# 将列表中的数值提取出来作为新列
merged_result_df <- merged_result_df %>%
  mutate(Chi_Square_P_Value_New = unlist(Chi_Square_P_Value)) %>%
  select(-Chi_Square_P_Value)

# 重命名新列
colnames(merged_result_df)[colnames(merged_result_df) == "Chi_Square_P_Value_New"] <- "Chi_Square_P_Value"

# 保存修改后的数据框到文件
write.csv(merged_result_df, "top0.01.dxy.HET_HOM.result.csv", row.names = FALSE)

