setwd("/fast3/group_crf/home/g20wangzhx36/m.hk/ismc/filter_DP_FROH/16inds");
options(scipen=999)

#load ROH segments:
arrHomFiles <- Sys.glob("/fast3/group_crf/home/g20wangzhx36/m.hk/ismc/filter_DP_FROH/16inds/ROH_segments/*.ROH.bed");
#remove M06 since it is same individual as Q1
#arrHomFiles <- arrHomFiles[!grepl("BTP_Q1", arrHomFiles)]
length(arrHomFiles)
datROHSegments <- NULL;
for (sF in arrHomFiles) {
  # 获取文件大小（以字节为单位）
  file_size <- file.info(sF)$size
  
  # 检查文件是否为空（文件大小为 0 字节表示为空）
  if (is.na(file_size) || file_size == 0) {
    next  # 如果文件为空，跳过本次循环，继续下一个文件
  }
  
  datThis <- read.table(sF, header = F, sep = "\t")
  datThis$Group <-  gsub('.ROH.bed', '', basename(sF), fixed = T)
  datROHSegments <- rbind(datROHSegments, datThis)
}

# 读取文件并指定列名
col_names <- c("chr", "homstart", "homend", "Group")
# 使用 names 函数指定列名
names(datROHSegments) <- col_names


unlink("rho.excludeROH.txt")

for (nChr in 1:23) {
  sF <- paste0(nChr, "/out.rho.10kb.bedgraph");
  if (!file.exists(sF)) {
    next;
  }
  dat20ind <- read.table(sF, header=T, sep="\t", fill=T, stringsAsFactors = F, check.names = FALSE)
  dat20ind$chromStart <- dat20ind$chromEnd - 9999;
  dat20ind$chromMid <- (dat20ind$chromStart + dat20ind$chromEnd)/2;
  dat19ind <- dat20ind[, colnames(dat20ind)!="BTP_Q0"] ; #EXCLUDE REF GENOME INDIVIDUAL

  dat19ind$chrom <- sChr <- paste0("mhkscf_",nChr);
  
  #set individual ROH segments to NA
  arrInd <- unique(datROHSegments$Group);
  arrInd <- arrInd[arrInd!="BTP_Q0"]
  for(sInd in arrInd  ) {
    datROHThis <- datROHSegments[datROHSegments$Group==sInd & datROHSegments$chr==sChr,]
    if (nrow(datROHThis) == 0) {
      next;
    }
    for(nROHRow in 1:nrow(datROHThis)) {
      dat19ind[dat19ind$chromMid >= datROHThis[nROHRow, 'homstart'] & dat19ind$chromMid <= datROHThis[nROHRow, 'homend'],  sInd ] <-NA;
    }
  }
  
  
  
  dat19ind$sample_mean <- rowMeans(dat19ind[,arrInd], na.rm = T)

  write.table(dat19ind[, c("chrom","chromStart" ,"chromEnd", "sample_mean")], file="rho.excludeROH.txt", sep="\t", quote = F, row.names = F, col.names = F,append = T)
  

}

#全基因组平均重组率0.09689947（已排除ROH区域）
