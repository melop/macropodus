setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp")

out <- "genetic_load_of_individual.Consurf9.allpops.refMHK.txt"

dat <- read.table("allpops.refMHK.NONSYN.polarized.out.filtered.withAD.dedup.withConsurf.dedup.bed", header=F, sep="\t")
dat <- dat[ , c(1:3, 6:ncol(dat) ) ]
#p_NA <- colSums(is.na(dat[, c(10:94)]))/nrow(dat)
#print(p_NA)
#dat <- dat[complete.cases(dat), ]
#dat_HIGH <- dat[dat$V10=="HIGH",]
dat_HIGH <- dat[dat$V190==9,]
p_NA <- colSums(is.na(dat_HIGH[, c(10:94)]))/nrow(dat_HIGH)
print(p_NA)

length_liftover <- 30740925 # liftover上的去除掉SV文件非双多态性的经过DP筛选的cds位点数 /data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/liftover/snp/filter_DP.R
#print(length_liftover)
#--------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------
# 自定义替换函数
replace_values <- function(x) {
  x <- ifelse(x == ".", "NA", x) # 如果是".", 则替换为"NA"
  
  # 将每个部分的"."替换为"0"
  x <- gsub("(?<=^|,)(\\.)(?=,|$)", "0", x, perl = TRUE)
  
  return(x)
}

# 提取子集并应用自定义函数
dat_HIGH[, 98:184] <- apply(dat_HIGH[, 98:184], 2, replace_values)
# 假设dat_HIGH是你的数据框
dat_HIGH <- dat_HIGH[!is.na(dat_HIGH[, ncol(dat_HIGH)]), ]
dat_HIGH <- dat_HIGH[!duplicated(dat_HIGH[, 1:3]), ]

# 定义一个函数，根据规则确定每个基因型对应的状态
assign_genotype_status <- function(genotype, reads_coverage, outgroup_AD) {
  if (is.na(genotype)) {
    return(NA)
  } else if (genotype == 0) {
    return("a")
  } else if (genotype == 1) {
    # 解析reads覆盖度
    reads_id <- strsplit(reads_coverage, ",")[[1]]
    if (length(reads_id) < 2) {
      # 处理reads_coverage缺失或格式错误的情况
      return(NA)
    }
    read_id1 <- as.numeric(reads_id[1])
    read_id2 <- as.numeric(reads_id[2])
    
    reads_out <- strsplit(outgroup_AD, ",")[[1]]
    if (length(reads_out) < 2) {
      # 处理outgroup_AD缺失或格式错误的情况
      return(NA)
    }
    read_out1 <- as.numeric(reads_out[1])
    read_out2 <- as.numeric(reads_out[2])
    
    if (is.na(read_id1) || is.na(read_id2) || is.na(read_out1) || is.na(read_out2)) {
      # 处理转换后为NA的情况
      return(NA)
    }
    
    if (read_out1 > read_out2){
      r <- "a"
      v <- "d"
    } else if (read_out1 < read_out2) {
      r <- "d"
      v <- "a"
    } else {
      r <- "equal"
      v <- "equal"
    }
    
    # 确定状态
    if (read_id1 > read_id2) {
      return(r)
    } else if (read_id1 < read_id2) {
      return(v)
    } else {
      # 如果相等，随机选择一个状态
      if (runif(1) > 0.5) {
        return(r)
      } else {
        return(v)
      }
    }
  } else if (genotype == 2) {
    return("d")
  }
}

# 初始化统计d的数量
d_counts <- rep(0, 85)

# 遍历每个位点的数据
for (i in 1:nrow(dat_HIGH)) {
  # 遍历每个个体的基因型和对应的reads覆盖度
  for (j in c(10:94)) {
    genotype <- dat_HIGH[i, j]
    reads_coverage <- dat_HIGH[i, j + 88]
    outgroup_AD <- dat_HIGH[i,184]
    
    # 确定状态
    status <- assign_genotype_status(genotype, reads_coverage, outgroup_AD)
    
    # 如果状态为d，增加对应的个体计数
    if (!is.na(status) && status == "d") {
      d_counts[j - 9] <- d_counts[j - 9] + 1
    }
    
    # 输出结果
    #cat(sprintf("Position %d, Sample %d: Genotype %s, Reads Coverage %s, Status %s\n",
    #            input_data[i, 1], j - 10, genotype, reads_coverage, status))
  }
}

# 统计每个个体所有位点随机抽取后 d 的数量
genetic_load <- d_counts
#genetic_load <- d_counts*2
print(genetic_load)

#个体遗传载荷标准化
length_filtered_SV <- nrow(dat) #筛选后位点数
print(length_filtered_SV)
#length_ref <- 429175390 #MHK参考基因组长度
genetic_load <- as.numeric(genetic_load)
standardized_genetic_load <- (genetic_load/(1-p_NA))/length_liftover
print(standardized_genetic_load)

# 保存结果为文本文件
result <- cbind(genetic_load, nrow(dat_HIGH), length_filtered_SV, p_NA, length_liftover, standardized_genetic_load)
colnames(result) <- c("genetic_load", "nrow_dat_HIGH", "length_filtered_SV", "percentage_NA","length_liftover", "standardized_genetic_load")
write.table(result, out, row.names = FALSE, quote = FALSE, sep = "\t")
