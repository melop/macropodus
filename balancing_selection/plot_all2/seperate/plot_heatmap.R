setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate")
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)

# 设定文件夹路径
path <- "/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate"

pop_list <- c("MHKhn", "MHKjx", "MHKmlh_qns", "MHKmls_pop1", "MHKmls_pop2", "MOPdsy", "MOPdxs_sg", "MOPfq_pt")

# 初始化存储结果的数据框
results <- data.frame(File = character(0),
                      Pop1 = character(0),
                      Pop2 = character(0),
                      Population1_Positive = numeric(0),
                      Population1_Negative = numeric(0),
                      Population2_Positive = numeric(0),
                      Population2_Negative = numeric(0),
                      P_Value = numeric(0))

for (i in 1:length(pop_list)) {
  for (j in i:length(pop_list)) {
    pop1 <- pop_list[i]
    pop2 <- pop_list[j]
    
    file1 <- file.path(path, paste(pop1, pop2, pop1, "seperated.bed", sep = "."))
    file2 <- file.path(path, paste(pop1, pop2, pop2, "seperated.bed", sep = "."))
    
    if (file.exists(file1) & file.exists(file2)) {
      data1 <- read.table(file1, header = FALSE, na.strings = "NA")
      data2 <- read.table(file2, header = FALSE, na.strings = "NA")
      
      colnames(data1) <- c("Chrom", "Start", "End", "NB")
      colnames(data2) <- c("Chrom", "Start", "End", "NB")
      
      result1 <- data1 %>%
        summarise(Population1_Positive = sum(NB > 0),
                  Population1_Negative = sum(NB == 0))
      
      result2 <- data2 %>%
        summarise(Population2_Positive = sum(NB > 0),
                  Population2_Negative = sum(NB == 0))
      
      # 卡方检验
      chisq_result <- chisq.test(matrix(c(result1$Population1_Positive, result1$Population1_Negative,
                                          result2$Population2_Positive, result2$Population2_Negative), ncol = 2))
      
      # 提取p值
      p_value <- chisq_result$p.value
      
      # 将结果保存到数据框
      results <- rbind(results, data.frame(File = paste(pop1, pop2, sep = "_"),
                                           Pop1 = pop1, Pop2 = pop2,
                                           Population1_Positive = result1$Population1_Positive,
                                           Population1_Negative = result1$Population1_Negative,
                                           Population2_Positive = result2$Population2_Positive,
                                           Population2_Negative = result2$Population2_Negative,
                                           P_Value = -log10(p_value)))
    }
  }
}

# 查看结果
print(results)
#write.table(results, file = "stats.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 创建一个空的矩阵，按照新顺序
mat <- matrix(NA, 
              nrow = length(pop_list),  # 新的 Pop1 数量
              ncol = length(pop_list),  # 新的 Pop2 数量
              dimnames = list(unique(results$Pop1), unique(results$Pop2)))  # 新的行列名

# 填充矩阵
for (i in 1:nrow(mat)) {
  for (j in 1:ncol(mat)) {
    pop1 <- unique(results$Pop1)[i]
    pop2 <- unique(results$Pop2)[j]
    value <- results$P_Value[results$Pop1 == pop1 & results$Pop2 == pop2]
    if (length(value) > 0) {
      mat[i, j] <- value
    } else {
      mat[i, j] <- NA
    }
  }
}

# 将矩阵转置，以使数据在对角线以下的位置也被填充
mat <- t(mat)

# 可选：将矩阵转换为数值型
mat <- as.matrix(mat)

# 打印矩阵
print(mat)

mark_mat <- ifelse(as.matrix(mat) >= 2, "**",  # P 值小于 0.01 标记为 **
                   ifelse(as.matrix(mat) >= 1.3010 & as.matrix(mat) < 2, "*", ""))  # P 值在 0.01 到 0.05 之间标记为 *


# 指定需要红色和蓝色显示的标签
red_labels <- c("MOPfq_pt", "MHKmlh_qns", "MHKjx")
blue_labels <- c("MHKhn", "MHKmls_pop1", "MHKmls_pop2", "MOPdsy", "MOPdxs_sg")

# 创建行标签的颜色映射
row_label_colors <- ifelse(rownames(mat) %in% red_labels, "red", 
                           ifelse(rownames(mat) %in% blue_labels, "blue", "black"))

# 创建列标签的颜色映射
col_label_colors <- ifelse(colnames(mat) %in% red_labels, "red", 
                           ifelse(colnames(mat) %in% blue_labels, "blue", "black"))

pdf(file = "heatmap_output2.pdf", width = 10, height = 8)  # 设置输出文件名和图像尺寸

# 创建热图
Heatmap(as.matrix(mat),
        cluster_rows = FALSE,  # 禁用行聚类
        cluster_columns = FALSE,  # 禁用列聚类
        col = colorRampPalette(c("white", "red"))(100),  # 设定颜色范围
        name = "-log10(P_Value)",  # 热图标题
        row_names_side = "left",  # 将行名放在左侧
        column_names_side = "top",  # 将列名放在顶部
        na_col = "#FFFFFF00",  # 设置 NA 的颜色为透明
        cell_fun = function(j, i, x, y, width, height, fill) {
          # 跳过 NA 单元格
          if (!is.na(as.matrix(mat)[i, j])) {
            grid.text(mark_mat[i, j], x, y)  # 在对应位置标记 * 或 **
          }
        },
        row_names_gp = gpar(col = row_label_colors, fontsize = 16),  # 设置行标签的颜色和字体大小
        column_names_gp = gpar(col = col_label_colors, fontsize = 16),  # 设置列标签的颜色和字体大小
        heatmap_legend_param = list(title = "-log10(p)", title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 16)),  # 设置图例标题和标签的字体大小
        row_title = "Population",  # 设置纵轴标题
        row_title_gp = gpar(fontsize = 20, fontface = "bold"),  # 设置纵轴标题的字体大小和样式
        column_title = "Population",  # 设置横轴标题
        column_title_gp = gpar(fontsize = 20, fontface = "bold")  # 设置横轴标题的字体大小和样式
)
# 关闭 PDF 设备
dev.off()
