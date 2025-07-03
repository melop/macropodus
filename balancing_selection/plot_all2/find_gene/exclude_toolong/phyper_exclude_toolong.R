#----------------------------------------------------------------------

# 设置工作目录
setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/find_gene/exclude_toolong")
library(data.table)
species_list <- c("MHKjx", "MHKmlh_qns", "MOPfq_pt")  # 请替换为实际的种群名称
N <- length(species_list)

process_species_pair <- function(sp1, sp2) {
  file_both_balanced <- paste0("../", sp1, ".", sp2, ".intersected.HET.balanced.bed")
  marker_file <- paste0(sp1, ".", sp2, ".BS_on_shearedHET.minDistance.txt")
  
  file_BS <- as.data.table(read.table(file_both_balanced, header = F, sep = "\t"))
  setnames(file_BS, 1:3, c("chr", "start", "end"))
  # 读取 marker 文件
  marker_data <- read.table(marker_file, header = TRUE, sep = "\t")
  # 筛选出长度不符合要求（小于 302 或大于 114959）的窗口
  invalid_markers <- marker_data[marker_data$minDistance < 302 | marker_data$minDistance > 114959, ]
  
  # 检查 invalid_markers 是否为空
  if (nrow(invalid_markers) == 0) {
    cat(paste0("组合 ", sp1, " 和 ", sp2, " 的 invalid_markers 为空，不进行筛选。\n"))
    valid_data_both_balanced <- file_BS
  } else {
    # 初始化一个空列表来存储分割后的结果
    split_list <- list()
    for (i in seq_along(invalid_markers$win)) {
      split_list[[i]] <- strsplit(as.character(invalid_markers$win[i]), split = "\\.")[[1]]
    }
  
    # 将列表转换为矩阵
    invalid_windows <- do.call(rbind, split_list)
    # 使用循环将第 2 列和第 3 列转换为 integer 类型
    for (col in 2:3) {
      invalid_windows[, col] <- as.numeric(invalid_windows[, col])
    }
    colnames(invalid_windows) <- c("chr", "start", "end")
    invalid_windows <- as.data.table(invalid_windows)
  
    # 创建一个唯一标识符列
    file_BS[, id := paste(chr, start, end, sep = "_")]
    invalid_windows[, id := paste(chr, start, end, sep = "_")]
  
    # 删除文件中与无效窗口相同的行
    valid_data_both_balanced <- file_BS[!id %in% invalid_windows$id]
    # 删除唯一标识符列
    valid_data_both_balanced[, id := NULL]
    }
  
  # 保存处理后的文件
  output_file <- paste0(sp1, ".", sp2, ".intersected.HET.balanced.exclude_toolong.bed")
  fwrite(valid_data_both_balanced, file = output_file, sep = "\t", col.names = FALSE)
}

for (i in 1:(N - 1)) {
  for (j in (i + 1):N) {
    sp1 <- species_list[i]
    sp2 <- species_list[j]
    cat(paste0("Processing combination: ", sp1, " and ", sp2, "\n"))
    process_species_pair(sp1, sp2)
  }
}


