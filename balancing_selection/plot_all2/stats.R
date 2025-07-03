setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2")
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
#library(gmodels)
#library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# 设定文件夹路径
path <- "/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2"

# 获取文件列表
files <- list.files(path, pattern = "\\.bothbeta\\.bed", full.names = TRUE)

# 初始化存储结果的数据框
results <- data.frame(File = character(0),
                      Pop1 = character(0),
                      Pop2 = character(0),
                      Population1_Positive = numeric(0),
                      Population1_Negative = numeric(0),
                      Population2_Positive = numeric(0),
                      Population2_Negative = numeric(0),
                      P_Value = numeric(0))

# 循环处理每个文件
for (file in files) {
  # 读取文件
  data <- read.table(file, header = FALSE, na.strings = "NA")
  
  # 提取种群和染色体编号
  file_name <- basename(file)
  pop <- gsub("\\.beta_on_intersectedHET\\.bed", "", file_name)
  pop <- strsplit(pop, "\\.")[[1]]
  pop1_id <- as.character(pop[1])
  pop2_id <- as.character(pop[2])
  
  # 更改列名
  colnames(data) <- c("BothHetInterval", "Pop1_nb", "Pop2_nb")
  
  # 统计各种群位点数大于0和等于0的行数
  result <- data %>%
    summarise(Population1_Positive = sum(Pop1_nb > 0),
              Population1_Negative = sum(Pop1_nb == 0),
              Population2_Positive = sum(Pop2_nb > 0),
              Population2_Negative = sum(Pop2_nb == 0)) 
  
  # 卡方检验
  chisq_result <- chisq.test(matrix(c(result$Population1_Positive, result$Population1_Negative, 
                                      result$Population2_Positive, result$Population2_Negative), ncol = 2))
  
  # 提取p值
  p_value <- chisq_result$p.value
  
  # 将结果保存到数据框
  results <- rbind(results, data.frame(Pop1 = pop1_id, Pop2 = pop2_id, result, P_Value = -log10(p_value)))
}

# 打印结果
print(results)
# 将结果保存到文件
#write.table(results, file = "stats.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#----------------------------------------------------------------------------------------------------------------
# 定义新的种群顺序
#new_order_Pop1 <- c("MHKhn", "MHKmls_pop1", "MHKmls_pop2","MOPdsy", "MOPdxs_sg", "MHKjx", "MHKmlh_qns", "MOPfq_pt")  # 新的 Pop1 排序
#new_order_Pop2 <- c("MHKhn", "MHKmls_pop1", "MHKmls_pop2","MOPdsy", "MOPdxs_sg", "MHKjx", "MHKmlh_qns", "MOPfq_pt")  # 新的 Pop2 排序

# 重新排列 Pop1 列
#results$Pop1 <- factor(results$Pop1, levels = new_order_Pop1)

# 重新排列 Pop2 列
#results$Pop2 <- factor(results$Pop2, levels = new_order_Pop2)

# 创建一个空的矩阵，按照新顺序
mat <- matrix(NA, 
              nrow = length(unique(results$Pop1)),  # 新的 Pop1 数量
              ncol = length(unique(results$Pop2)),  # 新的 Pop2 数量
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

pdf(file = "heatmap_output1.pdf", width = 10, height = 8)  # 设置输出文件名和图像尺寸

# 创建热图
Heatmap(as.matrix(mat),
        cluster_rows = FALSE,  # 禁用行聚类
        cluster_columns = FALSE,  # 禁用列聚类
        col = colorRampPalette(c("lightblue", "blue"))(100),  # 设定颜色范围
        name = "-log10(P_Value)",  # 热图标题
        row_names_side = "left",  # 将行名放在左侧
        column_names_side = "bottom",  # 将列名放在顶部
        na_col = "#FFFFFF00",  # 设置 NA 的颜色为透明
        cell_fun = function(j, i, x, y, width, height, fill) {
          # 跳过 NA 单元格
          if (!is.na(as.matrix(mat)[i, j])) {
            grid.text(mark_mat[i, j], x, y)  # 在对应位置标记 * 或 **
          }
        },
        row_names_gp = gpar(col = row_label_colors),  # 设置行标签的颜色
        column_names_gp = gpar(col = col_label_colors)  # 设置列标签的颜色
)
# 关闭 PDF 设备
dev.off()
