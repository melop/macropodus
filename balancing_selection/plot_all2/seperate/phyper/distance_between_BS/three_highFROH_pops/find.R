setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate/phyper/distance_between_BS/three_highFROH_pops")
# 读取三个文件
file_A <- read.table("MHKjx.BS.3pops_shearedHet.bed", header = FALSE)
file_B <- read.table("MHKmlh_qns.BS.3pops_shearedHet.bed", header = FALSE)
file_C <- read.table("MOPfq_pt.BS.3pops_shearedHet.bed", header = FALSE)

# 为列命名
colnames(file_A) <- c("chromosome", "start", "end")
colnames(file_B) <- c("chromosome", "start", "end")
colnames(file_C) <- c("chromosome", "start", "end")

# 定义一个函数来计算每个区间的中间位置
get_midpoint <- function(df) {
  df$midpoint <- (df$start + df$end) / 2
  return(df)
}

# 计算每个文件中区间的中间位置
file_A <- get_midpoint(file_A)
file_B <- get_midpoint(file_B)
file_C <- get_midpoint(file_C)

# 初始化各组合的最小距离和对应的位点信息
min_distance_AB <- Inf
closest_site_A_AB <- NULL
closest_site_B_AB <- NULL

min_distance_BC <- Inf
closest_site_B_BC <- NULL
closest_site_C_BC <- NULL

min_distance_AC <- Inf
closest_site_A_AC <- NULL
closest_site_C_AC <- NULL

# 计算 AB 之间的最短距离和位点
for (i in 1:nrow(file_A)) {
  site_A <- file_A[i, ]
  for (j in 1:nrow(file_B)) {
    site_B <- file_B[j, ]
    if (site_A$chromosome == site_B$chromosome) {
      distance_AB <- abs(site_A$midpoint - site_B$midpoint)
      if (distance_AB < min_distance_AB) {
        min_distance_AB <- distance_AB
        closest_site_A_AB <- site_A
        closest_site_B_AB <- site_B
      }
    }
  }
}

# 计算 BC 之间的最短距离和位点
for (i in 1:nrow(file_B)) {
  site_B <- file_B[i, ]
  for (j in 1:nrow(file_C)) {
    site_C <- file_C[j, ]
    if (site_B$chromosome == site_C$chromosome) {
      distance_BC <- abs(site_B$midpoint - site_C$midpoint)
      if (distance_BC < min_distance_BC) {
        min_distance_BC <- distance_BC
        closest_site_B_BC <- site_B
        closest_site_C_BC <- site_C
      }
    }
  }
}

# 计算 AC 之间的最短距离和位点
for (i in 1:nrow(file_A)) {
  site_A <- file_A[i, ]
  for (j in 1:nrow(file_C)) {
    site_C <- file_C[j, ]
    if (site_A$chromosome == site_C$chromosome) {
      distance_AC <- abs(site_A$midpoint - site_C$midpoint)
      if (distance_AC < min_distance_AC) {
        min_distance_AC <- distance_AC
        closest_site_A_AC <- site_A
        closest_site_C_AC <- site_C
      }
    }
  }
}

# 输出结果
cat("AB 之间最短距离:", min_distance_AB, "\n")
cat("文件 A 中对应的位点:\n")
print(closest_site_A_AB)
cat("文件 B 中对应的位点:\n")
print(closest_site_B_AB)

cat("\nBC 之间最短距离:", min_distance_BC, "\n")
cat("文件 B 中对应的位点:\n")
print(closest_site_B_BC)
cat("文件 C 中对应的位点:\n")
print(closest_site_C_BC)

cat("\nAC 之间最短距离:", min_distance_AC, "\n")
cat("文件 A 中对应的位点:\n")
print(closest_site_A_AC)
cat("文件 C 中对应的位点:\n")
print(closest_site_C_AC)
