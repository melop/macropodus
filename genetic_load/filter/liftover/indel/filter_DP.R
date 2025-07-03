setwd("/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/liftover/indel")
dat <- read.table("allpops.cds.bothDP.without_discarded.bed", header=F)
nrow(dat) #33832973
dat_filterDP <- dat[dat$V10>1244 & dat$V10<=5084 & dat$V14>1237 & dat$V14<=5079 & (dat$V10-dat$V14)>-43 & (dat$V10-dat$V14)<=43,]
nrow(dat_filterDP) #30741139
