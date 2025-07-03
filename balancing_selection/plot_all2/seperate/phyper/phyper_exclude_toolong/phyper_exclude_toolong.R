#----------------------------------------------------------------------

# 设置工作目录
setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/phyper/phyper_exclude_toolong")

# 加载必要的包
library(ComplexHeatmap)
library(circlize)
library(data.table)

# 假设 species_list 是包含所有种群名称的向量
species_list <- c("MHKhn", "MHKjx", "MHKmlh_qns", "MHKmls_pop1", "MHKmls_pop2", "MOPdsy", "MOPdxs_sg", "MOPfq_pt")  # 请替换为实际的种群名称
N <- length(species_list)

# 初始化 p 值矩阵
p_matrix <- matrix(NA, nrow = N, ncol = N)
rownames(p_matrix) <- species_list
colnames(p_matrix) <- species_list

# 定义一个函数来处理每一对种群组合
process_species_pair <- function(sp1, sp2) {
  # 构建文件路径
  file_both_balanced <- paste0("../", sp1, ".", sp2, ".", sp1, ".", sp2, ".bothbalanced.bed")
  file_only_sp1_balanced <- paste0("../", sp1, ".", sp2, ".", "only.", sp1, ".balanced.bed")
  file_only_sp2_balanced <- paste0("../", sp1, ".", sp2, ".", "only.", sp2, ".balanced.bed")
  file_both_no_balanced <- paste0("../", sp1, ".", sp2, ".", sp1, ".", sp2, ".bothno.balanced.bed")
  marker_file <- paste0("../distance_between_BS/", sp1, ".", sp2, ".BS_on_shearedHET.minDistance.txt")
  
  # 检查文件是否为空
  file_paths <- c(file_both_balanced, file_only_sp1_balanced, file_only_sp2_balanced, file_both_no_balanced, marker_file)
  empty_files <- sapply(file_paths, function(x) file.info(x)$size == 0)
  
  if (any(empty_files)) {
    if (file.info(marker_file)$size == 0) {
      cat(paste0("组合 ", sp1, " 和 ", sp2, " 的 marker 文件为空，筛选前后数据相同。\n"))
      # 统计筛选前的区间数量
      k_before <- nrow(read.table(file_both_balanced, header = FALSE, sep = "\t"))
      m_before <- k_before + nrow(read.table(file_only_sp1_balanced, header = FALSE, sep = "\t"))
      N_total_before <- k_before + nrow(read.table(file_only_sp1_balanced, header = FALSE, sep = "\t")) + 
        nrow(read.table(file_only_sp2_balanced, header = FALSE, sep = "\t")) + nrow(read.table(file_both_no_balanced, header = FALSE, sep = "\t"))
      n_before <- k_before + nrow(read.table(file_only_sp2_balanced, header = FALSE, sep = "\t"))
      
      # 计算筛选前的 p_value
      if (any(is.na(c(k_before, m_before, N_total_before, n_before)))) {
        p_value_before <- NA
      } else {
        p_value_before <- phyper(k_before, m_before, N_total_before - m_before, n_before, lower.tail = FALSE)
      }
      
      k <- k_before
      m <- m_before
      N_total <- N_total_before
      n <- n_before
      p_value <- p_value_before
    } else {
      cat(paste0("组合 ", sp1, " 和 ", sp2, " 的文件中存在空文件，跳过该组合。\n"))
      return(NULL)
    }
  } else {
    # 统计筛选前的区间数量
    k_before <- nrow(read.table(file_both_balanced, header = FALSE, sep = "\t"))
    m_before <- k_before + nrow(read.table(file_only_sp1_balanced, header = FALSE, sep = "\t"))
    N_total_before <- k_before + nrow(read.table(file_only_sp1_balanced, header = FALSE, sep = "\t")) + 
      nrow(read.table(file_only_sp2_balanced, header = FALSE, sep = "\t")) + nrow(read.table(file_both_no_balanced, header = FALSE, sep = "\t"))
    n_before <- k_before + nrow(read.table(file_only_sp2_balanced, header = FALSE, sep = "\t"))
    
    # 计算筛选前的 p_value
    if (any(is.na(c(k_before, m_before, N_total_before, n_before)))) {
      p_value_before <- NA
    } else {
      p_value_before <- phyper(k_before, m_before, N_total_before - m_before, n_before, lower.tail = FALSE)
    }
    
    # 检查 marker 文件是否存在
    if (!file.exists(marker_file)) {
      cat(paste0("在处理组合 ", sp1, " 和 ", sp2, " 时，marker 文件 ", marker_file, " 不存在，行数视为 0。\n"))
      k <- 0
      m <- 0
      N_total <- 0
      n <- 0
    } else {
      # 读取 marker 文件
      marker_data <- read.table(marker_file, header = TRUE, sep = "\t")
      # 筛选出长度不符合要求（小于 302 或大于 114959）的窗口
      invalid_markers <- marker_data[marker_data$minDistance < 302 | marker_data$minDistance > 114959, ]
      
      # 检查 invalid_markers 是否为空
      if (nrow(invalid_markers) == 0) {
        cat(paste0("组合 ", sp1, " 和 ", sp2, " 的 invalid_markers 为空，筛选前后数据相同。\n"))
        k <- k_before
        m <- m_before
        N_total <- N_total_before
        n <- n_before
        p_value <- p_value_before
      } else {
        # 初始化一个空列表来存储分割后的结果
        split_list <- list()
        for (i in seq_along(invalid_markers$win)) {
          split_list[[i]] <- strsplit(as.character(invalid_markers$win[i]), split = "\\.")[[1]]
        }
        
        # 将列表转换为矩阵
        invalid_windows <- do.call(rbind, split_list)
        
        # 检查 invalid_windows 的行数和列数
        if (!is.null(dim(invalid_windows)) && nrow(invalid_windows) > 0 && ncol(invalid_windows) >= 3) {
          # 使用循环将第 2 列和第 3 列转换为 integer 类型
          for (col in 2:3) {
            invalid_windows[, col] <- as.numeric(invalid_windows[, col])
          }
          colnames(invalid_windows) <- c("chr", "start", "end")
          invalid_windows <- as.data.table(invalid_windows)
        } else {
          cat(paste0("在处理组合 ", sp1, " 和 ", sp2, " 时，invalid_windows 行数或列数不足，不进行后续类型转换操作。\n"))
          invalid_windows <- NULL  # 将 invalid_windows 设置为 NULL
        }
        
        # 读取要统计的文件
        file_data_both_balanced <- as.data.table(read.table(file_both_balanced, header = FALSE, sep = "\t"))
        file_data_only_sp1_balanced <- as.data.table(read.table(file_only_sp1_balanced, header = FALSE, sep = "\t"))
        file_data_only_sp2_balanced <- as.data.table(read.table(file_only_sp2_balanced, header = FALSE, sep = "\t"))
        file_data_both_no_balanced <- as.data.table(read.table(file_both_no_balanced, header = FALSE, sep = "\t"))
        
        setnames(file_data_both_balanced, 1:3, c("chr", "start", "end"))
        setnames(file_data_only_sp1_balanced, 1:3, c("chr", "start", "end"))
        setnames(file_data_only_sp2_balanced, 1:3, c("chr", "start", "end"))
        setnames(file_data_both_no_balanced, 1:3, c("chr", "start", "end"))
        
        if (!is.null(invalid_windows) && nrow(invalid_windows) > 0) {
          # 创建一个唯一标识符列
          file_data_both_balanced[, id := paste(chr, start, end, sep = "_")]
          invalid_windows[, id := paste(chr, start, end, sep = "_")]
          
          # 删除文件中与无效窗口相同的行
          valid_data_both_balanced <- file_data_both_balanced[!id %in% invalid_windows$id]
          k <- nrow(valid_data_both_balanced)
          
          file_data_only_sp1_balanced[, id := paste(chr, start, end, sep = "_")]
          valid_data_only_sp1_balanced <- file_data_only_sp1_balanced[!id %in% invalid_windows$id]
          
          file_data_only_sp2_balanced[, id := paste(chr, start, end, sep = "_")]
          valid_data_only_sp2_balanced <- file_data_only_sp2_balanced[!id %in% invalid_windows$id]
          
          file_data_both_no_balanced[, id := paste(chr, start, end, sep = "_")]
          valid_data_both_no_balanced <- file_data_both_no_balanced[!id %in% invalid_windows$id]
          
          # 删除唯一标识符列
          valid_data_both_balanced[, id := NULL]
          valid_data_only_sp1_balanced[, id := NULL]
          valid_data_only_sp2_balanced[, id := NULL]
          valid_data_both_no_balanced[, id := NULL]
          
          m <- k + nrow(valid_data_only_sp1_balanced)
          N_total <- k + nrow(valid_data_only_sp1_balanced) + nrow(valid_data_only_sp2_balanced) + nrow(valid_data_both_no_balanced)
          n <- k + nrow(valid_data_only_sp2_balanced)
        } else {
          k <- nrow(file_data_both_balanced)
          m <- k + nrow(file_data_only_sp1_balanced)
          N_total <- k + nrow(file_data_only_sp1_balanced) + nrow(file_data_only_sp2_balanced) + nrow(file_data_both_no_balanced)
          n <- k + nrow(file_data_only_sp2_balanced)
        }
        
        # 进行 phyper 检验
        if (any(is.na(c(k, m, N_total, n)))) {
          p_value <- NA
        } else {
          p_value <- phyper(k, m, N_total - m, n, lower.tail = FALSE)
        }
      }
    }
  }
  
  # 存储筛选前和筛选后的信息
  result <- list(
    before = list(k = k_before, m = m_before, N_total = N_total_before, n = n_before, p_value = p_value_before),
    after = list(k = k, m = m, N_total = N_total, n = n, p_value = p_value)
  )
  
  return(result)
}

# 遍历所有种群组合
all_results <- list()
for (i in 1:(N - 1)) {
  for (j in (i + 1):N) {
    sp1 <- species_list[i]
    sp2 <- species_list[j]
    cat(paste0("Processing combination: ", sp1, " and ", sp2, "\n"))
    result <- process_species_pair(sp1, sp2)
    if (!is.null(result)) {
      all_results[[paste(sp1, sp2, sep = "_")]] <- result
      p_value <- result$after$p_value
      p_matrix[i, j] <- p_value
      p_matrix[j, i] <- p_value
    }
  }
}

# 输出所有结果
print(all_results)
# 输出 p 值矩阵
print(p_matrix)
# 将矩阵转换为数据框（如果矩阵有行名和列名，数据框也会保留）
mat_df <- as.data.frame(p_matrix, stringsAsFactors = FALSE)
#-----------------------------------------------------------------------------绘图
pop_list <- c("MHKhn", "MHKjx", "MHKmlh_qns", "MHKmls_pop1", "MHKmls_pop2", "MOPdsy", "MOPdxs_sg", "MOPfq_pt")

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
pdf(file = "Balance_map_Het.phyper.exclude_toolong.pdf", width = 10, height = 8)

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
        na_col = "gray",  # 设置 NA 的颜色为透明
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

