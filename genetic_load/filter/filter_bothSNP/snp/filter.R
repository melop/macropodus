setwd("/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_bothSNP/snp")
dat <- read.table("allpops.cds.bothDP.bothAF.bed", header=F, sep="\t")
nrow(dat) 
dat <- dat[dat$V8>1244 & dat$V8<=5084 & dat$V10>1237 & dat$V10<=5079 & (dat$V8-dat$V10)>-43 & (dat$V8-dat$V10)<=43,]
nrow(dat) 
dat$Fixed1 <- (dat$V14==0 | dat$V14==1); #单个数值用||，数组用|
dat$Fixed2 <- (dat$V18==0 | dat$V18==1);

# 筛选出Fixed1和Fixed2相等的行
dat_bothSNP <- dat[(dat$Fixed1 == dat$Fixed2) & !is.na(dat$Fixed1) & !is.na(dat$Fixed2), ]
nrow(dat_bothSNP)

# 计算dat的补集，包括NA行
dat_non_bothSNP <- dat[!(dat$Fixed1 == dat$Fixed2 & !is.na(dat$Fixed1) & !is.na(dat$Fixed2)), ]
nrow(dat_non_bothSNP)

dat1 <- dat_bothSNP[,c(1:3)]
dat2 <- dat_bothSNP[,c(4:6)]
dat3 <- dat_non_bothSNP[,c(1:3)]
dat4 <- dat_non_bothSNP[,c(4:6)]

write.table(dat1, file = "allpops.refMHK.NONSYN.filtered.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dat2, file = "allpops.refMOP.NONSYN.filtered.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dat3, file = "allpops.refMHK.NONSYN.discarded.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dat4, file = "allpops.refMOP.NONSYN.discarded.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

