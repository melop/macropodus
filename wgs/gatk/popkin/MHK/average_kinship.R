setwd("/data2/projects/zwang/m.hk/Structure/Popkin")
dat <- read.csv("MHK.kinship.tsv", header = T, sep = "\t")

pops <- read.csv("MHK.sample_names.csv", header = F, sep = ",")

# 检查文件B第二列的长度是否和文件A的行数、列数一致
if (nrow(dat) != nrow(pops) || ncol(dat) != nrow(pops)) {
  stop("文件B第二列的长度与文件A的行数和列数不一致，无法设置行名和列名。")
}

# 用文件B的第二列设置文件A的行名和列名
# 生成唯一的行名
unique_names <- make.unique(as.character(pops[, 2]))

# 检查 unique_names 的长度是否和 dat 的行数一致
if (nrow(dat) != length(unique_names)) {
  stop("生成的唯一行名数量与矩阵的行数不一致。")
}

# 设置行名
rownames(dat) <- unique_names
colnames(dat) <- pops[, 2]


# 获取所有不同的种群名称
populations <- unique(colnames(dat))

# 初始化一个空向量，用于存储每个种群的平均亲缘系数
average_kinship <- numeric(length(populations))
names(average_kinship) <- populations


# 遍历每个种群
for (pop in populations) {
  # 找出列名属于当前种群的索引
  col_indices <- which(colnames(dat) == pop)
  
  # 找出将行名去掉后缀后属于当前种群的索引
  row_indices <- which(sub("\\.\\d+$", "", rownames(dat)) == pop)
  
  # 如果该种群只有一个个体，则平均亲缘系数设为 NA
  if (length(row_indices) <= 1 || length(col_indices) <= 1) {
    average_kinship[pop] <- NA
    next
  }
  
  # 提取该种群个体之间的亲缘系数矩阵的上三角部分（不包括对角线）
  sub_matrix <- dat[row_indices, col_indices]
  upper_tri <- sub_matrix[upper.tri(sub_matrix)]
  
  # 计算平均亲缘系数
  average_kinship[pop] <- mean(upper_tri)
}

# 输出结果
average_kinship

# 提取整个矩阵的上三角部分（不包括对角线）
all_upper_tri <- dat[upper.tri(dat)]

# 计算所有种群的平均亲缘系数（包括种群间个体的亲缘系数）
if (length(all_upper_tri) > 0) {
  overall_average_kinship <- mean(all_upper_tri)
  print("所有种群的平均亲缘系数（包括种群间）：")
  print(overall_average_kinship)
} else {
  print("没有有效的亲缘系数数据用于计算所有种群的平均亲缘系数。")
  overall_average_kinship <- NA
}

# 将结果保存到文件
output_df <- data.frame(
  Population = c(names(average_kinship), "Overall"),
  Average_Kinship = c(average_kinship, overall_average_kinship)
)

write.table(output_df, "MHK.average_kinship.txt", sep = "\t", na = "NA", quote = FALSE, row.names = FALSE)




#-----------------------------------------------------------------------之前用的/data2/projects/dyao/macropodus/mhk/kinship/MHK/kinship.tsv
# 计算所有个体的平均亲缘系数
average_kinship <- mean(dat[lower.tri(dat)])
print(average_kinship)


# 判断 dat 的类型，如果不是矩阵，则尝试转换为矩阵
if (!is.matrix(dat)) {
  dat <- as.matrix(dat)
}

# 获取对角线上的值
diagonal_values <- diag(dat)

# 计算对角线上值的平均值
average_diagonal <- mean(diagonal_values)
print(average_diagonal)

#MHKpopkin = 0.9633936
#MOPpopkin = 0.8854923