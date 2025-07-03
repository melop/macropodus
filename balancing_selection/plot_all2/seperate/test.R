setwd("/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate")

dat <- read.csv("MHKjx.MHKmlh_qns.intersected.HET.bed", header = F, sep = "\t")

dat$V4 <- dat$V3-dat$V2
