setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/indel/ROH")

#--------------------------------------------------------------------------------------
# 文件列表，假设所有需要处理的文件都在一个列表中
file_list <- list.files(path = "/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/indel/ROH", pattern = "allpops\\.refMHK\\.SV\\.filtered\\.(.*)\\.(ROH|HET)\\.bed$", full.names = TRUE)

# 提取个体名
individuals <- unique(gsub("allpops\\.refMHK\\.SV\\.filtered\\.(.*?)\\.(ROH|HET)\\.bed$", "\\1", basename(file_list)))

# 计算区间长度的函数
calculate_interval_length <- function(dat) {
  return(sum(dat$V3 - dat$V2))
}
# 计算遗传载荷的函数
calculate_genetic_load <- function(dat) {
  dat <- dat[ , c(1:3, 6:ncol(dat) ) ]
  dat_HIGH <- dat[dat$V10 == "HIGH", ]
  p_NA <- colSums(is.na(dat_HIGH[, c(10:ncol(dat_HIGH))]))/nrow(dat_HIGH)
  print(p_NA)
  
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
  
  random_dat_HIGH <- matrix(NA, nrow = nrow(dat_HIGH), ncol = ncol(dat_HIGH) - 9)
  for (i in 10:ncol(dat_HIGH)) {
    for (j in 1:nrow(dat_HIGH)) {
      random_dat_HIGH[j, i - 9] <- random_selection(dat_HIGH[j, i])
    }
  }
  
  genetic_load <- colSums(random_dat_HIGH == "d", na.rm = TRUE)
  return(genetic_load/(1-p_NA))
}


# 遍历每个个体
for (individual in individuals) {
  # 构建ROH和HET文件名
  roh_file <- paste0("allpops.refMHK.SV.filtered.", individual, ".ROH.bed")
  het_file <- paste0("allpops.refMHK.SV.filtered.", individual, ".HET.bed")
  
  # 完整路径
  roh_path <- file.path("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/indel/ROH", roh_file)
  het_path <- file.path("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/indel/ROH", het_file)
  
  # 区间长度文件路径
  interval_roh_file <- paste0(individual, ".ROH.bed") # 根据实际情况修改
  interval_het_file <- paste0(individual, ".HET.bed") # 根据实际情况修改
  
  # 完整路径
  interval_roh_path <- file.path("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/ROH", interval_roh_file)
  interval_het_path <- file.path("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/ROH", interval_het_file)
  
  # 读取区间长度
  interval_roh_dat <- read.table(interval_roh_path, header = FALSE, sep = "\t")
  interval_het_dat <- read.table(interval_het_path, header = FALSE, sep = "\t")
  length_roh <- calculate_interval_length(interval_roh_dat)
  length_het <- calculate_interval_length(interval_het_dat)
  
  # 读取ROH和HET数据
  if (file.exists(roh_path) && file.exists(het_path)) {
    roh_dat <- read.table(roh_path, header = FALSE, sep = "\t")
    het_dat <- read.table(het_path, header = FALSE, sep = "\t")
    
    # 计算遗传载荷
    genetic_load_roh <- calculate_genetic_load(roh_dat)
    genetic_load_het <- calculate_genetic_load(het_dat)
    
    # 计算遗传载荷比值（先除以区间长度）
    normalized_genetic_load_roh <- genetic_load_roh / length_roh
    normalized_genetic_load_het <- genetic_load_het / length_het
    # 将标准化的遗传载荷保存为第一列和第二列，并包含列名
    output_data <- data.frame(
      normalized_genetic_load_roh = normalized_genetic_load_roh,
      normalized_genetic_load_het = normalized_genetic_load_het
    )
    output_file <- paste0("genetic_load_", individual, ".txt")
    write.table(output_data, file = output_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    
    # 打印结果
    print(paste("Processed:", individual))
    print(output_data)
  } else {
    print(paste("Files for individual", individual, "are missing."))
  }
}
#--------------------------------------------------------------------------------------

