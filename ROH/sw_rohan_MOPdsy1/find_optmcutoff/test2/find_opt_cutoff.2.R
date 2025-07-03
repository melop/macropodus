# 设置工作目录，根据实际情况调整路径
setwd("/data2/projects/zwang/m.hk/ROH/sw_rohan_DSY1/find_optmcutoff/test2")

# 定义不同测序深度文件的文件名前缀列表
depth_prefixes <- c("X10", "X12", "X14", "X16", "X18", "X2", "X20", "X22", "X24", "X4", "X6", "X8")

# 用于存储每个文件对应的结果（准确率等信息）
all_results <- data.frame(Depth = character(), Accuracy = numeric(), TP = numeric(), TN = numeric(), FP = numeric(), FN = numeric(), estroh = numeric(),stringsAsFactors = FALSE)

datROHTruth <- read.table("roh.bed", sep="\t", header=F);
datROHTruth$len <- datROHTruth$V3 - datROHTruth$V2

# 固定的阈值
nCutoff <- 0.00082

# 循环处理每个测序深度的文件
for(depth_prefix in depth_prefixes) {
  # 构建文件名
  het_file_name <- paste0(depth_prefix, ".rohan.het.bed")
  win_file_name <- paste0(depth_prefix, ".roh.wins.bed")
  
  # 读取杂合位点文件
  datROHAN <- read.table(het_file_name, sep = "\t", header = F, stringsAsFactors = F)
  datROHAN <- datROHAN[complete.cases(datROHAN), ]
  datROHAN <- datROHAN[datROHAN$V1 %in% unique(datROHTruth$V1) , ];
  
  colnames(datROHAN) <- c('chr', 'winstart', 'winend', 'esthet')
  
  # 读取纯合区域交集文件
  datROHWin <- read.table(win_file_name, sep = "\t", header = F)
  datROHWin <- datROHWin[complete.cases(datROHWin),]
  
  # 按照固定阈值计算相关数值
  nROH_on_ROH <- 0;
  nROH_on_HET <- 0;
  nHET_on_ROH <- 0;
  nHET_on_HET <- 0;
  
  
  #  for(sChr in arrChrs) {
  datROHTrueOnChr <- datROHTruth; #  datROHTruth[datROHTruth$V1 == sChr, ];
  datROHANonChr <- datROHAN; # datROHAN[datROHAN$chr == sChr, ];
  datOvlpOnChr <- datROHWin; # datROHWin[datROHWin$V1 == sChr, ];
  
  datROHTrueOnChr$coord <- paste(datROHTrueOnChr$V1, datROHTrueOnChr$V2, datROHTrueOnChr$V3, sep=":");
  datROHANonChr$coord <- paste(datROHANonChr$chr, datROHANonChr$winstart, datROHANonChr$winend, sep=":");
  datOvlpOnChr$coord <- paste(datOvlpOnChr$V1, datOvlpOnChr$V2, datOvlpOnChr$V3 ,sep=":");
  
  datCalledROH <- datROHANonChr[datROHANonChr$esthet<=nCutoff, ];
  
  datCalledNonROH <- datROHANonChr[datROHANonChr$esthet > nCutoff, ];
  
  datCalledROHOvlp <- datOvlpOnChr[datOvlpOnChr$V4<=nCutoff, ];
  
  datCalledNonROHOvlp <- datOvlpOnChr[datOvlpOnChr$V4>nCutoff, ];
  
  nROHOvlpLen <- sum(datCalledROHOvlp$V8);
  nROHTrueLen <- sum(datROHTrueOnChr$len);
  nChrLen <- max(datROHANonChr$winend);
  
  nROH_on_ROH <- nROH_on_ROH + nROHOvlpLen;
  
  datROHonHet <- datCalledROH[!(datCalledROH$coord %in% datCalledROHOvlp$coord), ]
  nROHonHetLenOnChr <- 0;
  if (nrow(datROHonHet) > 0) {
    nROHonHetLenOnChr <- sum(datROHonHet$winend - datROHonHet$winstart + 1);
  }
  nROHonHetLenOnChr <- nROHonHetLenOnChr + sum(datCalledROHOvlp$V3-datCalledROHOvlp$V2+1 - datCalledROHOvlp$V8)
  nROH_on_HET <- nROH_on_HET + nROHonHetLenOnChr
  
  nHetOnROHChr <- sum(datCalledNonROHOvlp$V8);
  nHET_on_ROH <- nHET_on_ROH + nHetOnROHChr;
  
  nHetOnHetChr <- sum(datCalledNonROHOvlp$V3-datCalledNonROHOvlp$V2+1-datCalledNonROHOvlp$V8);
  datHetOnHet <- datCalledNonROH[!(datCalledNonROH$coord %in% datCalledROHOvlp$coord), ]
  nHET_on_HET <- nHET_on_HET + nHetOnHetChr + sum(datHetOnHet$winend - datHetOnHet$winstart + 1)
  
  nEstFROH <- (nROH_on_ROH + nROH_on_HET) / (nROH_on_ROH + nROH_on_HET + nHET_on_ROH + nHET_on_HET);
  
  correct <- (nROH_on_ROH + nHET_on_HET) / (nROH_on_ROH + nROH_on_HET + nHET_on_ROH + nHET_on_HET)
  
  TP <-  nROH_on_ROH  / (nROH_on_ROH + nHET_on_ROH)
  TN <-  nHET_on_HET / (nROH_on_HET + nHET_on_HET)
  FP <-  nROH_on_HET/ (nROH_on_HET + nHET_on_HET)
  FN <-  nHET_on_ROH/ (nROH_on_ROH + nHET_on_ROH)
  
  # 将当前文件的结果添加到总结果数据框中
  all_results <- rbind(all_results, data.frame(Depth = depth_prefix, Accuracy =  correct, TP= TP, TN= TN, FP= FP, FN= FN, estroh= nEstFROH))
}

print(all_results)

# 指定保存结果的txt文件名，可根据需求修改
output_file_name <- "fixed_threshold_results.1e-3.9.txt"
# 使用write.table函数将结果数据框保存到txt文件中，sep="\t"表示用制表符分隔列，row.names = FALSE表示不保存行名
write.table(all_results, file = output_file_name, sep = "\t", row.names = FALSE)

library(ggplot2)
# 从文件名前缀中提取测序深度数值（假设前缀中的数字就是测序深度值）
extract_depth <- function(prefix) {
  depth_num <- as.numeric(gsub("X", "", prefix))
  return(depth_num)
}
# 提取测序深度数值到新列（假设前缀中的数字就是测序深度值）
all_results$DP <- sapply(all_results$Depth, extract_depth)

# 使用ggplot2绘制曲线图
p <- ggplot(all_results, aes(x = DP, y = Accuracy)) +
  geom_line() +
  labs(x = "Dpeth", y = "Accuracy", title = "Accuracy of cutoff 0.00082 in different depth") +
  theme_bw()+
  theme(
    # 设置标题字体大小，单位为pt（磅），这里设置为16pt，可根据需求调整
    plot.title = element_text(size = 20),
    # 设置坐标轴标签字体大小，这里设置为14pt，可根据需求调整
    axis.title = element_text(size = 20),
    # 设置坐标轴刻度标签字体大小，这里设置为12pt，可根据需求调整
    axis.text = element_text(size = 20),
    # 设置图形与边框的边距，单位为pt（磅），依次为上、右、下、左边距，可根据需求调整数值
    plot.margin = unit(c(20, 20, 20, 20), "pt")
  )

# 指定保存的pdf文件名及路径，可根据实际需求修改
pdf_file_name <- "DSY-1.fixed_shreshold.accuracy_curve.1e-3.9.pdf"
# 使用ggsave函数将绘制好的图保存为pdf文件，width和height参数可调整图片尺寸
ggsave(pdf_file_name, plot = p, width = 10, height = 6)
