setwd("/fast3/group_crf/home/g20wangzhx36/m.hk/Structure/Pixy/output/0.99")

dxy <- read.table("../pixy_output_dxy.txt", header = T, sep = "\t")
fst <- read.table("../pixy_output_fst.txt", header = T, sep = "\t")

dxy <- dxy[as.numeric(substr(dxy$chromosome, 8, nchar(dxy$chromosome))) < 24, ]
fst <- fst[as.numeric(substr(fst$chromosome, 8, nchar(fst$chromosome))) < 24, ]


# 先筛选出一个种群为MHK开头，另一个为MOP开头的行
sub_dxy <- dxy[(grepl("^MHK", dxy$pop1) & grepl("^MOP", dxy$pop2)) |
                 (grepl("^MOP", dxy$pop1) & grepl("^MHK", dxy$pop2)), ]

sub_fst <- fst[(grepl("^MHK", fst$pop1) & grepl("^MOP", fst$pop2)) |
                 (grepl("^MOP", fst$pop1) & grepl("^MHK", fst$pop2)), ]

# 先筛除第6列为NA的行
sub_dxy <- sub_dxy[!is.na(sub_dxy$avg_dxy), ]
sub_fst <- sub_fst[!is.na(sub_fst$avg_wc_fst), ]


# 创建一个空列表来存储不同种群组合的筛选结果
filtered_results_dxy <- list()
filtered_results_fst <- list()

# 获取所有不同的种群组合（在已筛选出的MHK和MOP组合数据中）
unique_combinations <- unique(paste(sub_dxy$pop1, sub_dxy$pop2, sep = "-"))

#---------------------dxy提取前1%，fst提取后1%---------------------------------

# 遍历每个种群组合
for (combination in unique_combinations) {
  pop1 <- strsplit(combination, "-")[[1]][1]
  pop2 <- strsplit(combination,"-")[[1]][2]
#  print(pop1)
#  print(pop2)
  subset_dxy <- sub_dxy[(sub_dxy$pop1 == pop1 &
                        sub_dxy$pop2 == pop2), ]

  
  # 计算当前组合中需要筛选出的前1%的行数
  top_dxy_rows <- ceiling(nrow(subset_dxy) * 0.01)
  
  # 按照dxy值对当前组合的数据进行降序排序
  sorted_subset_dxy <- subset_dxy[order(-subset_dxy$avg_dxy), ]
  
  # 筛选出dxy值在前1%的行
  filtered_dxy <- sorted_subset_dxy[1:top_dxy_rows, ]
  
  # 将当前组合的筛选结果添加到列表中
  filtered_results_dxy[[combination]] <- filtered_dxy
  
  # 根据当前种群组合筛选数据
  subset_fst <- sub_fst[(sub_fst$pop1 == pop1 &
                           sub_fst$pop2 == pop2), ]
  
  # 计算当前组合中需要筛选出的前1%的行数
  top_fst_rows <- ceiling(nrow(subset_fst) * 0.01)
  
  # 按照fst值对当前组合的数据进行升序排序
  sorted_subset_fst <- subset_fst[order(subset_fst$avg_wc_fst), ]
  
  # 筛选出fst值在后1%的行
  filtered_fst <- sorted_subset_fst[1:top_fst_rows, ]
  
  # 将当前组合的筛选结果添加到列表中
  filtered_results_fst[[combination]] <- filtered_fst
}

# 将列表中的结果合并成一个数据框（如果需要的话）
final_dxy <- do.call(rbind, filtered_results_dxy)
final_fst <- do.call(rbind, filtered_results_fst)

test_unique_combinations <- unique(paste(final_dxy$pop1, final_dxy$pop2, sep = "-"))

#-------------------------------------------------------------------------------------共1110种窗口。42种组合中dxy均在各自前1%的窗口有87个10kb窗口
# 根据第三、四、五列创建一个新的列，用于标识每条记录的组合信息
final_dxy$window_identifier <- paste(final_dxy$chromosome, final_dxy$window_pos_1, final_dxy$window_pos_2, sep = "_")

window_count_dxy <- table(final_dxy$window_identifier)

common_windows_dxy <- final_dxy[final_dxy$window_identifier %in% names(window_count_dxy[window_count_dxy == 42]), ] 

write.table(common_windows_dxy, file = "pairedMHKMOP.dxy.top0.01.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#--------------------------------------------------------------------------------------共1301种窗口。42种组合中fst均在各自后1%的窗口有2个10kb窗口
# 根据第三、四、五列创建一个新的列，用于标识每条记录的组合信息
final_fst$window_identifier_fst <- paste(final_fst$chromosome, final_fst$window_pos_1, final_fst$window_pos_2, sep = "_")

window_count_fst <- table(final_fst$window_identifier)

common_windows_fst <- final_fst[final_fst$window_identifier %in% names(window_count_fst[window_count_fst == 42]), ]

write.table(common_windows_fst, file = "pairedMHKMOP.fst.bottom.0.01.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#--------------------------------------------------------------------------------

# 提取MHKmlh_qns和MOPfq_pt组合的数据
subset_fst_MOPfq_pt_MHKmlh_qns <- final_fst[(final_fst$pop1 == "MHKmlh_qns" & final_fst$pop2 == "MOPfq_pt") |
                                              (final_fst$pop1 == "MOPfq_pt" & final_fst$pop2 == "MHKmlh_qns"), ]

# 提取MHKjx和MOPfq_pt组合的数据
subset_fst_MOPfq_pt_MHKjx <- final_fst[(final_fst$pop1 == "MHKjx" & final_fst$pop2 == "MOPfq_pt") |
                                         (final_fst$pop1 == "MOPfq_pt" & final_fst$pop2 == "MHKjx"), ]
# 检查是否有重合的窗口
overlapping_windows <- intersect(subset_fst_MOPfq_pt_MHKmlh_qns$window_identifier_fst, subset_fst_MOPfq_pt_MHKjx$window_identifier_fst)
print(overlapping_windows)
