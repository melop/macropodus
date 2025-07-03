setwd("/data2/projects/zwang/m.hk/ROH/sw_rohan_DSY1")
dat <- read.table("var.stats.mhkscf_1.tsv", header=T, sep="\t");
hist(as.numeric(dat$DP))

# 计算DP分布的0.01和0.99分位数
dp_quantiles <- quantile(as.numeric(dat$DP), probs = c(0.01, 0.99), na.rm = TRUE)

# 输出结果
print(dp_quantiles)
#1% 6
#99% 457