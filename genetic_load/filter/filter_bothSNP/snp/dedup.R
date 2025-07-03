setwd("/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_bothSNP/snp")

ori_dat <- read.table("allpops.refMHK.NONSYN.withAF.bed", header = F, sep = "\t")

unique_data <- unique(ori_dat)

# 筛选重复的区间
duplicates <- ori_dat[duplicated(ori_dat) | duplicated(ori_dat, fromLast = TRUE), ] #12546 sites

write.table(unique_data, "allpops.cds.depth.dedup.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
