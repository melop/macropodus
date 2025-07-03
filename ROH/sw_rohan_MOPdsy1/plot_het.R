#setwd("D:\\raycui\\projects\\Macropodus_hongkongensis\\gatk\\slidingwin_het");
setwd("/data2/projects/zwang/m.hk/ROH/sw_rohan_DSY1")
sOut <- "mhk_on_mop.hom.win.tsv"
arrChr <- paste0("mopscf_" , 1:23);
datHet <- read.table( gzfile('DSY-1.50000.100.het.tsv.gz','rt') , header=F);
# 设置全局字体大小
par(cex = 1.5)  # 整体缩放因子，可根据需要调整该值
# 绘制 log10(datHet$V9) 的直方图
hist(log10(datHet$V9), breaks = 1000,
     main = "MOPdsy-1",
     xlab = "Log10 (heterozygosity)",
     ylab = "Frequency",
     cex.main = 1.5,
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlim = c(-5.0, -1.0),
)
hist(datHet$V9, breaks=100000,xlim=c(0,0.001))
hist(log10(datHet$V9), breaks=1000)
# 在 x = -3.9 处添加一条红色垂直线
abline(v = -3.9, col = "red", lwd = 3)
#hist_result <- hist(datHet$V9, breaks = 100000, plot = FALSE)
# 计算 log10(datHet$V9)
log_values <- log10(datHet$V9)

# 统计小于 -3.9 和大于 3.9 的数据个数
less_than_neg_3_9 <- sum(log_values < -3.9)
greater_than_neg_3_9 <- sum(log_values > -3.9)

# 计算数据总数
total_count <- length(log_values)

# 计算比例
proportion_less_than_neg_3_9 <- less_than_neg_3_9 / total_count
proportion_greater_than_neg_3_9 <- greater_than_neg_3_9 / total_count

# 输出结果
cat("小于 -3.9 的数据所占比例：", proportion_less_than_neg_3_9, "\n")
cat("大于 3.9 的数据所占比例：", proportion_greater_than_3_9, "\n")




# 找到频数最高的V9值
#most_frequent_v9 <- hist_result$mids[which.max(hist_result$counts)]

#most_frequent_v9


nHomCutoff <- 3.162278e-05 #if at this cutoff, count as homozygous there is no significant difference between refMHK and refMOP, so nHomCutoff and nHomCutoff2 use refMHK's.
nHomCutoff2 <- 5e-5 ; #if at this cutoff, count as ambiguous
nMinEffSites <- 10000; #count as ambiguous if lower than this number of effective sites
nMaxAmbgStretch <- 1000000; #max ambiguous stretch before breaking block.

datHomWin <- NULL;
write.table(data.frame(chr='chr', ambstart='ambstart', homstart='homstart',  homend='homend', ambend='ambend', stringsAsFactors = F), file = sOut , row.names = F, col.names = F, sep="\t", quote = F);
for(sChr in arrChr) {
  datHetChr <- datHet[datHet$V1 == sChr, ];

  nPrevState <- -1;
  nAmbStart <- -1;
  nAmbEnd <- -1;
  nHomStart <- -1;
  nHomEnd <- -1;
  nPrevEnd <- -1;
  nPrevStart <- -1;
  for(i in 1:nrow(datHetChr)) {
    nPosLeft <- datHetChr[i, 2];
    nPosRight <-datHetChr[i, 3];
    nEff <- datHetChr[i, 5];
    nTheta <- datHetChr[i, 9];
    nState <- 1; #het 
    if (nEff < nMinEffSites || (nTheta <=nHomCutoff2 && nTheta > nHomCutoff) ) {
      nState <- 2;
    } else if (nTheta <= nHomCutoff) {
      nState <- 0;
    }
    
    bBreakAmb <- F;

    if (nPrevState == -1) {
      nPrevState <- nState;
      nPrevStart <- 1;
      nPrevEnd <- 1;
    }
    
    if (nState ==2) {
      if (nPrevState == nState) {
        if ((nPosRight - nAmbStart+1) > nMaxAmbgStretch) {
          bBreakAmb <- T;
          if (nHomStart ==-1) {
            nAmbStart <- nPosLeft;
          }
        }
      } else if (nPrevState == 1) {
        nAmbStart <- nPosLeft;
      } else if (nPrevState == 0) {
        nHomEnd <- nPrevEnd;
        nAmbEnd <- nPosRight;
      }
      
    }
    
    if (nState == 0) {
      if (nPrevState == nState) {
        nHomEnd <- nPosRight;
      } else if (nPrevState == 2) {
        if (nAmbStart < nHomStart) {
          nHomEnd <- nPosRight;
        } else {
          nHomStart <- nPosLeft;
          nHomEnd <- nPosRight;
        }
      } else if (nPrevState == 1) {
       # cat("0 to 1:", nPosLeft, "\n");
        nAmbStart <- nPosLeft;
        nHomStart <- nPosLeft;
      }
    }
    
    if (nState == 1) {
      if (nPrevState == nState) {
        
      } else if (nPrevState == 2) {
        nAmbEnd <- nPrevEnd;
        bBreakAmb <- T;
      } else if (nPrevState == 0) {
        nAmbEnd <- nPrevEnd;
        nHomEnd <- nPrevEnd;
        bBreakAmb <- T;
      }
    }
    
    if (bBreakAmb) {
      if (nHomStart > 0 && nAmbEnd >= nHomEnd) {
        #datHomWin <- rbind(datHomWin, data.frame(chr=sChr, ambstart=nAmbStart, homstart=nHomStart,  homend=nHomEnd, ambend=nAmbEnd, stringsAsFactors = F));
        if (nAmbStart == -1) {
          nAmbStart <- nHomStart;
        }
        write.table(data.frame(chr=sChr, ambstart=nAmbStart, homstart=nHomStart,  homend=nHomEnd, ambend=nAmbEnd, stringsAsFactors = F), file = sOut , append = T, row.names = F, col.names = F, sep="\t", quote = F);
        
      }
      if (nAmbStart != nPosLeft) {
      nAmbStart <- -1;
      }
      nAmbEnd <- -1;
      nHomStart <- -1;
      nHomEnd <- -1;
      
    }
    nPrevState <- nState;
    nPrevEnd <- nPosRight;
    nPrevStart <- nPosLeft;
    
  }
}
