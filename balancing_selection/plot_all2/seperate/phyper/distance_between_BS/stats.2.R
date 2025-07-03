# 设置工作目录
setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/phyper/distance_between_BS")

# 定义处理文件对的函数
process_file_pair <- function(file_A_path, file_B_path) {
  # 读取文件 A 和文件 B
  file_A <- read.table(file_A_path, header = FALSE, col.names = c("chromosome", "start", "end", "win_chr", "win_start", "win_end"))
  file_B <- read.table(file_B_path, header = FALSE, col.names = c("chromosome", "start", "end", "win_chr", "win_start", "win_end"))
  
  # 按染色体和窗口分组
  groups_A <- split(file_A, list(file_A$chromosome, file_A$win_start, file_A$win_end))
  groups_B <- split(file_B, list(file_B$chromosome, file_B$win_start, file_B$win_end))
  
  # 初始化一个列表来存储每个窗口的最短距离
  window_min_distances <- list()
  
  # 遍历所有可能的染色体 - 窗口组合
  for (key in intersect(names(groups_A), names(groups_B))) {
    current_A <- groups_A[[key]]
    current_B <- groups_B[[key]]
    
    # 检查当前分组是否有数据
    if (nrow(current_A) > 0 && nrow(current_B) > 0) {
      # 初始化当前窗口的最短距离为无穷大
      window_min_distance <- Inf
      
      # 计算当前染色体 - 窗口组合中 A 和 B 位点之间的距离
      for (i in 1:nrow(current_A)) {
        for (j in 1:nrow(current_B)) {
          distance <- abs(current_A$start[i] - current_B$start[j])
          
          # 如果当前距离小于当前窗口的最短距离，更新最短距离
          if (distance < window_min_distance) {
            window_min_distance <- distance
          }
        }
      }
      
      # 将当前窗口的最短距离存入列表
      window_min_distances[[key]] <- window_min_distance
      
      # 打印当前窗口的最短距离
      cat("文件对:", basename(file_A_path), "和", basename(file_B_path), "的窗口", key, "的最短距离:", window_min_distance, "\n")
    }
  }
  
  return(window_min_distances)
}


# 定义文件列表
file_list <- c("MHKhn", "MHKjx", "MHKmlh_qns", "MHKmls_pop1", "MHKmls_pop2", "MOPdsy", "MOPdxs_sg", "MOPfq_pt")

# 遍历所有可能的文件对组合
for (i in seq_along(file_list)) {
  for (j in i:length(file_list)) {
    file_A_name <- paste0(file_list[i], ".BS.", file_list[i], ".", file_list[j], ".shearedHET.bed")
    file_B_name <- paste0(file_list[j], ".BS.", file_list[i], ".", file_list[j], ".shearedHET.bed")
    
    # 检查文件是否存在
    if (file.exists(file_A_name) && file.exists(file_B_name)) {
      # 处理当前文件对
      result <- process_file_pair(file_A_name, file_B_name)
      
      # 生成输出文件名
      output_file_name <- paste0(file_list[i], ".", file_list[j], ".BS_on_shearedHET.minDistance.txt")
      
      # 将结果转换为数据框
      result_df <- data.frame(
        win = names(result),
        minDistance = unlist(result)
      )
      
      # 将数据框写入文件
      write.table(result_df, file = output_file_name, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
      cat("已将文件对", file_A_name, "和", file_B_name, "的结果保存到", output_file_name, "\n")
    } else {
      cat("文件对", file_A_name, "和", file_B_name, "中存在文件不存在，跳过处理。\n")
    }
  }
}

