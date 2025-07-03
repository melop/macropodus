setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/indel")

out <- "genetic_load_of_individual.allpops.refMHK.txt"

#dat_origin <- read.table("/data2/projects/zwang/test/DeepVatiant/MOPhn/test/filter_bothSNP_for_individual/indel_on_liftover/MOPhn.refMHK.SV.w.indelmodifier.polarized.out.bed", header=F, sep="\t")
#p_NA <- colSums(is.na(dat_origin[, c(12:ncol(dat_origin))]))/nrow(dat_origin)
#print(p_NA)

dat <- read.table("/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_bothSNP/indel/allpops.refMHK.SV.polarized.out.filtered.bed", header=F, sep="\t")
dat <- dat[ , c(1:3, 6:ncol(dat) ) ]
#p_NA <- colSums(is.na(dat[, c(10:ncol(dat))]))/nrow(dat)
#print(p_NA)
#dat <- dat[complete.cases(dat), ]
dat_HIGH <- dat[dat$V10=="HIGH",]
p_NA <- colSums(is.na(dat_HIGH[, c(10:ncol(dat_HIGH))]))/nrow(dat_HIGH)
print(p_NA)

length_liftover <- 30741139 # liftover上的去除掉SV文件非双多态性的经过DP筛选的cds位点数 /data2/projects/zwang/test/DeepVatiant/DPcutoff_for_all/liftover_for_individual/indel_forGATK/filter_DP.R
#print(length_liftover)

# 定义函数，根据基因型从相应的等位基因中随机抽取一个
random_selection <- function(genotype) {
  if (is.na(genotype)) {
    return(NA)
  } else if (genotype == "0") {
    return("a")
  } else if (genotype == "1") {
    return(sample(c("a", "d"), 1))
  } else if (genotype == "2") {
    return("d")
  }
}
# 随机抽取每个个体每个位点的等位基因
random_dat_HIGH <- matrix(NA, nrow = nrow(dat_HIGH), ncol = ncol(dat_HIGH)-9)
for (i in 10:ncol(dat_HIGH)) {
  for (j in 1:nrow(dat_HIGH)) {
    random_dat_HIGH[j, i-9] <- random_selection(dat_HIGH[j, i])
  }
}

# 统计每个个体所有位点随机抽取后 d 的数量
genetic_load <- colSums(random_dat_HIGH == "d", na.rm = TRUE)
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

