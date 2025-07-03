setwd("/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_DP/refMOP")

ori_dat <- read.table("allpops.cds.depth.bed", header = F, sep = "\t")

unique_data <- unique(ori_dat)

# 检查第三列减去第二列的值是否不等于1的行
# 计算差值
diffs <- unique_data$V3 - unique_data$V2

# 找到差值不等于1的行
invalid_rows <- unique_data[diffs != 1, ]

# 筛选重复的区间
duplicates <- ori_dat[duplicated(ori_dat) | duplicated(ori_dat, fromLast = TRUE), ] #10781 sites


write.table(unique_data, "allpops.cds.depth.dedup.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
