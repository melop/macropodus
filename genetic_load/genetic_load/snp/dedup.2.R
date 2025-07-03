setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp")

ori_dat <- read.table("allpops.refMHK.NONSYN.polarized.out.filtered.withAD.dedup.withConsurf.bed", header = F, sep = "\t")

unique_data <- ori_dat[!duplicated(ori_dat[, 1:3]) & !duplicated(ori_dat[, 1:3], fromLast = TRUE), ]

# 筛选重复的区间
duplicates <- ori_dat[duplicated(ori_dat[, 1:3]) | duplicated(ori_dat[, 1:3], fromLast = TRUE), ] #154 sites

write.table(unique_data, "allpops.refMHK.NONSYN.polarized.out.filtered.withAD.dedup.withConsurf.dedup.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

