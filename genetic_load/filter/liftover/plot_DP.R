setwd("/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/liftover")

my_data <- read.table("allpops.cds.bothDP.bed", header = F)
dat <- my_data[!duplicated(my_data[,1:3]),]

##DP_refMHK
# 计算2.5%和97.5%的置信区间
confidence_interval <- quantile(dat$V10, c(0.025, 0.975))

hist(dat$V10, main = "DP of cds", xlab = "DP", ylab = "Frequency", breaks = "FD")

# 添加标记
abline(v = confidence_interval[1], col = "red", lty = 2, lwd = 2)
abline(v = confidence_interval[2], col = "red", lty = 2, lwd = 2)
text(confidence_interval[1], 100000, paste("2.5% CI:", round(confidence_interval[1], 2)), adj = c(1, 0.5), col = "red")
text(confidence_interval[2], 100000, paste("97.5% CI:", round(confidence_interval[2], 2)), adj = c(0, 0.5), col = "red")
#deduplicated dat: 2.5%~1244 97.5%~5084


##DP_refMOP
# 计算2.5%和97.5%的置信区间
confidence_interval <- quantile(dat$V14, c(0.025, 0.975))

hist(dat$V14, main = "DP of cds", xlab = "DP", ylab = "Frequency", breaks = "FD")

# 添加标记
abline(v = confidence_interval[1], col = "red", lty = 2, lwd = 2)
abline(v = confidence_interval[2], col = "red", lty = 2, lwd = 2)
text(confidence_interval[1], 100000, paste("2.5% CI:", round(confidence_interval[1], 2)), adj = c(1, 0.5), col = "red")
text(confidence_interval[2], 100000, paste("97.5% CI:", round(confidence_interval[2], 2)), adj = c(0, 0.5), col = "red")
#deduplicated dat: 2.5%~1237 97.5%~5079

##DP_refMHK-DP_refMOP
#difference <- abs(dat[, 10] - dat[, 17])
difference <- dat[, 10] - dat[, 14]


hist(difference, xlab="DPvariation", col="lightblue", main="Histogram of DPvariation", breaks = "FD")

confidence_interval <- quantile(difference,c(0.025,0.975))

abline(v = confidence_interval[1], col = "red", lty = 2, lwd = 2)
abline(v = confidence_interval[2], col = "red", lty = 2, lwd = 2)
text(confidence_interval[1], 1000000, paste("2.5% CI:", round(confidence_interval[1], 2)), adj = c(1, 0.5), col = "red")
text(confidence_interval[2], 1000000, paste("97.5% CI:", round(confidence_interval[2], 2)), adj = c(0, 0.5), col = "red")

#deduplicated dat: 2.5%~-43 97.5%~43
