setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/Ration_of_HomHet/HIGH")

dat <- read.table("//data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/allpops.refMHK.NONSYN.polarized.out.filtered.withAD.dedup.withConsurf.dedup.bed", header=F, sep="\t")
dat <- dat[ , c(1:3, 6:ncol(dat) ) ]
#p_NA <- colSums(is.na(dat[, c(10:ncol(dat))]))/nrow(dat)
#print(p_NA)
#dat <- dat[complete.cases(dat), ]
dat_HIGH <- dat[dat$V10=="HIGH",]

# 第10列到最后一列是个体的基因型数据
genotype_data <- dat_HIGH[, 10:94]

# 将genotype_data中的数据转化为数值类型
genotype_data <- as.data.frame(lapply(genotype_data, function(x) as.numeric(as.character(x))))

# 统计每个个体中基因型为1的位点数
count_1 <- apply(genotype_data, 2, function(x) sum(x == 1, na.rm = TRUE))

# 统计每个个体中基因型为2的位点数
count_2 <- apply(genotype_data, 2, function(x) sum(x == 2, na.rm = TRUE))

count_1_2 <- count_1 + count_2

# 计算每个个体中基因型为2的占比
proportion <- count_2 / count_1_2


# 输出结果
proportion

# 将 proportion_1 和 proportion_2 合并为一个数据框
result <- data.frame(
  Individual = colnames(genotype_data),  # 每个个体的名称
  Proportion = proportion           # 基因型1和2的比例
)

# 将结果保存为 txt 文件
write.table(result, file = "genotype_proportions.HIGH.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#---------------------------------------------------------------------------------
#检查是否存在重复位点
# 创建唯一标识符，把 start 和 end 组合成字符串
#dat_HIGH$position <- paste(dat_HIGH$V1, dat_HIGH$V2, sep = "-")

# 统计每个位置的出现次数
#position_counts <- table(dat_HIGH$position)

# 统计重复出现的位点数量
#repeated_positions <- sum(position_counts > 1)

# 统计所有位点的总数
#total_positions <- length(position_counts)

# 计算重复位点所占的比例
#proportion_repeated <- repeated_positions / total_positions

# 输出结果
#cat("重复出现的位点在所有位点中所占的比例为：", proportion_repeated, "\n")

#---------------------------------------------------------------------------------
