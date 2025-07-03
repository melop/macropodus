setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/ROH/HIGH")

#--------------------------------------------------------------------------------------
# 文件列表，假设所有需要处理的文件都在一个列表中
file_list <- list.files(path = "/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/ROH/HIGH", pattern = "allpops\\.refMHK\\.NONSYN\\.filtered\\.(.*)\\.(ROH|HET)\\.bed$", full.names = TRUE)

# 提取个体名
individuals <- unique(gsub("allpops\\.refMHK\\.NONSYN\\.filtered\\.(.*?)\\.(ROH|HET)\\.bed$", "\\1", basename(file_list)))

# 计算区间长度的函数
calculate_interval_length <- function(dat) {
  return(sum(dat$V3 - dat$V2))
}

# 自定义替换函数
replace_values <- function(x) {
  x <- ifelse(x == ".", "NA", x) # 如果是".", 则替换为"NA"
  
  # 将每个部分的"."替换为"0"
  x <- gsub("(?<=^|,)(\\.)(?=,|$)", "0", x, perl = TRUE)
  
  return(x)
}

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

# 计算遗传载荷的函数
calculate_genetic_load <- function(dat) {
  dat <- dat[ , c(1:3, 6:ncol(dat) ) ]
  dat_HIGH <- dat[dat$V10=="HIGH", ]
  p_NA <- colSums(is.na(dat_HIGH[, c(10:94)]))/nrow(dat_HIGH)
  print(p_NA)
  
  # 提取子集并应用自定义函数
  dat_HIGH[, 98:184] <- apply(dat_HIGH[, 98:184], 2, replace_values)
  # 假设dat_HIGH是你的数据框
  dat_HIGH <- dat_HIGH[!is.na(dat_HIGH[, ncol(dat_HIGH)]), ]
  dat_HIGH <- dat_HIGH[!duplicated(dat_HIGH[, 1:3]), ]
  
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
    }
  }
  genetic_load <- d_counts
  return(genetic_load/(1-p_NA))
}


# 遍历每个个体
for (individual in individuals) {
  # 构建ROH和HET文件名
  roh_file <- paste0("allpops.refMHK.NONSYN.filtered.", individual, ".ROH.bed")
  het_file <- paste0("allpops.refMHK.NONSYN.filtered.", individual, ".HET.bed")
  
  # 完整路径
  roh_path <- file.path("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/ROH/HIGH", roh_file)
  het_path <- file.path("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/ROH/HIGH", het_file)
  
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

