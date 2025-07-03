setwd("/fast3/group_crf/home/g20wangzhx36/m.hk/Structure/Pixy/output")

dat <- read.table("dxy_average.txt", header = F, sep = " ")

# 提取所有涉及的种群名字
all_populations <- unique(c(dat[, 1], dat[, 2]))

# 创建一个空矩阵，大小根据种群数量确定
mat <- matrix(nrow = length(all_populations),
              ncol = length(all_populations),
              dimnames = list(all_populations, all_populations))


# 填充矩阵，先按照原始顺序放置值
for (i in 1:nrow(dat)) {
  row_name <- dat[i, 1]
  col_name <- dat[i, 2]
  mat[row_name, col_name] <- dat[i, 5]
}

# 将所有值统一放置到pop1vspop2对应的位置
for (i in 1:length(all_populations)) {
  for (j in 1:length(all_populations)) {
    if (i!= j) {
      pop1 <- all_populations[i]
      pop2 <- all_populations[j]
      value <- mat[pop2, pop1]
      if (!is.na(value)) {
        mat[pop1, pop2] <- value
        mat[pop2, pop1] <- NA
      }
    }
  }
}

# 此时mat就是转换好的矩阵，你可以根据需要进行后续操作，比如查看矩阵内容
print(mat)
# 将处理好的矩阵保存为CSV文件
write.csv(mat, file = "paired.averageDxy.csv", row.names = TRUE)
