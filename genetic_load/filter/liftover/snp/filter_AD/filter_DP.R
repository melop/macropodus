setwd("/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/liftover/snp/filter_AD")
dat <- read.table("allpops.cds.bothDP.without_discarded.bed", header=F)
nrow(dat) #35055454
dat_filterDP <- dat[dat$V10>1220 & dat$V10<=5103 & dat$V14>1209 & dat$V14<=5094 & (dat$V10-dat$V14)>-60 & (dat$V10-dat$V14)<=62,]
nrow(dat_filterDP) #31936383
