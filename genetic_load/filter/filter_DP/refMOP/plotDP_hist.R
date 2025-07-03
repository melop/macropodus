setwd("/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_DP/refMOP")

my_data <- read.table("allpops.cds.depth.bed", header = F)

##DP_refMHK
# 计算2.5%和97.5%的置信区间
confidence_interval <- quantile(my_data$V4, c(0.025, 0.975))

hist(my_data$V4, main = "DP of cds", xlab = "DP", ylab = "Frequency", breaks = "FD")

# 添加标记
abline(v = confidence_interval[1], col = "red", lty = 2, lwd = 2)
abline(v = confidence_interval[2], col = "red", lty = 2, lwd = 2)
text(confidence_interval[1], 100000, paste("2.5% CI:", round(confidence_interval[1], 2)), adj = c(1, 0.5), col = "red")
text(confidence_interval[2], 100000, paste("97.5% CI:", round(confidence_interval[2], 2)), adj = c(0, 0.5), col = "red")
#2.5%～，97.5%～
