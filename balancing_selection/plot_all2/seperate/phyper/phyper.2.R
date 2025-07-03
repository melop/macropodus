setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/phyper")

# 加载ComplexHeatmap包
library(ComplexHeatmap)
library(circlize)

count_intervals <- function(file_path) {
  # 使用readLines函数读取文件内容为字符向量
  lines <- readLines(file_path)
  # 获取字符向量的长度，即文件的行数，也就是杂合区间数量
  count <- length(lines)
  return(count)
}


# 假设list是包含文件名的向量
list <- c("MHKhn", "MHKjx", "MHKmlh_qns", "MHKmls_pop1", "MHKmls_pop2", "MOPdsy", "MOPdxs_sg", "MOPfq_pt")
N <- length(list)
p_matrix <- matrix(NA, nrow = N, ncol = N)
# 使用种群名设置矩阵的行名和列名
rownames(p_matrix) <- list
colnames(p_matrix) <- list
for (i in 1:(N - 1)) {
  for (j in (i + 1):N) {
    file_both_balanced <- paste(list[i], list[j], list[i], list[j], "bothbalanced.bed", sep = ".")
    file_only_sp1_balanced <- paste(list[i], list[j], "only", list[i], "balanced.bed", sep = ".")
    file_only_sp2_balanced <- paste(list[i], list[j], "only", list[j], "balanced.bed", sep = ".")
    file_both_no_balanced <- paste(list[i], list[j], list[i], list[j], "bothno.balanced.bed", sep = ".")
    
    # 统计区间数量
    k <- count_intervals(file_both_balanced)
    m <- count_intervals(file_both_balanced)+count_intervals(file_only_sp1_balanced)
    N_total <- count_intervals(file_both_balanced)+count_intervals(file_only_sp1_balanced)+count_intervals(file_only_sp2_balanced)+count_intervals(file_both_no_balanced)
    n <- count_intervals(file_both_balanced)+count_intervals(file_only_sp2_balanced)
    
    # 进行phyper检验
    p_value <- phyper(k, m, N_total - m, n, lower.tail = FALSE)
#    print(paste("For combination", list[i], "and", list[j], "p - value:", p_value))
    p_matrix[i, j] <- p_value
    p_matrix[j, i] <- p_value
  }
}
print(p_matrix)

# 将矩阵转换为数据框（如果矩阵有行名和列名，数据框也会保留）
mat_df <- as.data.frame(p_matrix, stringsAsFactors = FALSE)

# 使用write.csv函数保存为csv文件，设置quote = FALSE表示字符串不用双引号括起来
#write.csv(mat_df, file = "Balance_map_Het.phyper_test.csv", quote = FALSE)

#------------------------------------------------------------------------------
# 定义包含种群名的向量（与矩阵的行名和列名对应）
pop_list <- c("MHKhn", "MHKjx", "MHKmlh_qns", "MHKmls_pop1", "MHKmls_pop2", "MOPdsy", "MOPdxs_sg", "MOPfq_pt")

# 定义原始矩阵（这里直接按照你提供的数据赋值，实际应用中可能需要从其他地方读取等方式获取）
#p_matrix <- matrix(c(NA, 4.582563e-10, 5.675116e-16, 0.033598835, 7.000106e-06, 2.089918e-05, 1.238213e-05, 9.797527e-07,
#                     4.582563e-10, NA, 1.165695e-05, 0.003497897, 1.897776e-06, 6.471781e-08, 3.069368e-14, 2.177649e-02,
#                     5.675116e-16, 1.165695e-05, NA, 0.236066607, 4.599058e-18, 2.554032e-12, 9.905040e-06, 2.358868e-04,
#                     0.033598835, 0.003497897, 0.236066607, NA, 0.056848895, 0.232109718, 0.001302579, 0.350212267,
#                     7.000106e-06, 1.897776e-06, 4.599058e-18, 0.056848895, NA, 2.589352e-06, 7.739010e-04, 2.647496e-05,
#                     2.089918e-05, 6.471781e-08, 2.554032e-12, 0.232109718, 2.589352e-06, NA, 2.601626e-07, 2.942600e-07,
#                     1.238213e-05, 3.069368e-14, 9.905040e-06, 0.001302579, 7.739010e-04, 2.601626e-07, NA, 4.123712e-07,
#                     9.797527e-07, 2.177649e-02, 2.358868e-04, 0.350212267, 2.647496e-05, 2.942600e-07, 4.123712e-07, NA),
#                   nrow = 8, byrow = TRUE,
#                   dimnames = list(pop_list, pop_list))

# 只保留矩阵的左下半部分（不包含对角线，可根据需求调整是否包含对角线）的数据
new_p_matrix <- matrix(NA, nrow = nrow(p_matrix), ncol = nrow(p_matrix), dimnames = dimnames(p_matrix))
for (i in 1:nrow(p_matrix)) {
  for (j in 1:i) {
    new_p_matrix[i, j] = p_matrix[i, j]
  }
}

mat <- as.matrix(new_p_matrix)
# 指定需要红色和蓝色显示的标签
red_labels <- c("MOPfq_pt", "MHKmlh_qns", "MHKjx")
blue_labels <- c("MHKhn", "MHKmls_pop1", "MHKmls_pop2", "MOPdsy", "MOPdxs_sg")

# 创建行标签的颜色映射
row_label_colors <- ifelse(rownames(mat) %in% red_labels, "red", 
                           ifelse(rownames(mat) %in% blue_labels, "blue", "black"))

# 创建列标签的颜色映射
col_label_colors <- ifelse(colnames(mat) %in% red_labels, "red", 
                           ifelse(colnames(mat) %in% blue_labels, "blue", "black"))


# 创建pdf文件用于保存热图，设置输出文件名和图像尺寸
pdf(file = "Balance_map_Het.phyper_test.2.pdf", width = 10, height = 8)

# 对矩阵中的值进行 -log10转换（处理0值情况，添加极小值避免取对数出现 -Inf）
mat_transformed <- -log10(mat + 1e-100)
# 获取 -log10转换后矩阵中的最小值和最大值，用于设定颜色映射范围
min_value <- min(mat_transformed, na.rm = TRUE)
max_value <- max(mat_transformed, na.rm = TRUE)
# 根据最小值和最大值创建自定义的颜色映射函数，确保颜色渐变更贴合数据范围
col_fun <- colorRamp2(c(min_value, max_value), c("white", "red"))

# 创建热图
Heatmap(mat_transformed,
        cluster_rows = FALSE,  # 禁用行聚类
        cluster_columns = FALSE,  # 禁用列聚类
        col = col_fun,  # 使用自定义的颜色映射函数，实现根据 -log10值从白到红渐变
        name = "-log10(P_Value)",  # 热图标题
        row_names_side = "left",  # 将行名放在左侧
        column_names_side = "top",  # 将列名放在顶部
        na_col = "#FFFFFF00",  # 设置 NA 的颜色为透明
        cell_fun = function(j, i, x, y, width, height, fill) {
          # 判断p值范围并进行相应标记
          if (!is.na(mat[i, j]) && mat[i, j] < 0.01) {
            grid.text("**", x, y)
          } else if (!is.na(mat[i, j]) && mat[i, j] < 0.05) {
            grid.text("*", x, y)
          }
        },
        row_names_gp = gpar(col = row_label_colors, fontsize = 24),  # 设置行标签的颜色和字体大小
        column_names_gp = gpar(col = col_label_colors, fontsize = 24),  # 设置列标签的颜色和字体大小
        heatmap_legend_param = list(title = "-log10(p)", title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 24)),  # 设置图例标题和标签的字体大小
        row_title = "Population",  # 设置纵轴标题
        row_title_gp = gpar(fontsize = 26, fontface = "bold"),  # 设置纵轴标题的字体大小和样式
        column_title = "Population",  # 设置横轴标题
        column_title_gp = gpar(fontsize = 26, fontface = "bold")  # 设置横轴标题的字体大小和样式
)
dev.off()  # 关闭图形设备，确保pdf文件正确保存