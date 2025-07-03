#setwd("D:\\raycui\\projects\\Macropodus_hongkongensis\\gatk\\slidingwin_het");
#sOut <- "mhk.hom.win.ver2.tsv"
#sOut <- "mhk.hom.win.ver3.tsv"
#sOut <- "mhk.hom.win.ver4.tsv"
#sOut <- "mhk.hom.win.ver5.tsv"
#sOut <- "mhk.hom.win.ver6.tsv"
sOut <- "mhk.hom.win.ver7.tsv"
arrChr <- paste0("mhkscf_" , 1:23);
datHet <- read.table( gzfile('mhk.50000.100.het.tsv.gz','rt') , header=F);
#hist(log10(datHet$V9), breaks=1000)
#nHomCutoff <- 3.162278e-05 #if at this cutoff, count as homozygous
#nHomCutoff <- 5e-05;
#nHomCutoff <- 5e-04;
#nHomCutoff <- 8e-04;
nHomCutoff <- 1e-03;
#nHomCutoff <- 0.000158489 # 1e-3.8

datHomWin <- NULL;
write.table(data.frame(chr='chr', ambstart='ambstart', homstart='homstart',  homend='homend', ambend='ambend', stringsAsFactors = F), file = sOut , row.names = F, col.names = F, sep="\t", quote = F);
for(sChr in arrChr) {
  datHetChr <- datHet[datHet$V1 == sChr, ];

  nPrevState <- -1;
  nHomStart <- -1;
  nHomEnd <- -1;
  nPrevMid <- -1;
  bBreakAmb <- F;
  for(i in 1:nrow(datHetChr)) {
    nPosMid <- datHetChr[i, 4];
#    nEff <- datHetChr[i, 5];
    nTheta <- datHetChr[i, 9];
    nState <- 1; #het 
    if (nTheta <= nHomCutoff) {
      nState <- 0;
    }
    

    if (nPrevState == -1) {
      nPrevState <- nState;
      nPrevMid <- 1;
    }
    
    if (nState == 0) {
      if (nPrevState == nState) {
	nHomEnd <- nPosMid;
      } else if (nPrevState == 1) {
       # cat("0 to 1:", nPosLeft, "\n");
        nHomStart <- nHomEnd <- nPosMid;
      }
    }
    
    if (nState == 1) {
      if (nPrevState == nState) {
      } else if (nPrevState == 0) {
        bBreakAmb <- T;
      }
    }
    
    if (bBreakAmb) {
      if (nHomStart > 0 ) {
        write.table(data.frame(chr=sChr, ambstart=nHomStart, homstart=nHomStart,  homend=nHomEnd, ambend=nHomEnd, stringsAsFactors = F), file = sOut , append = T, row.names = F, col.names = F, sep="\t", quote = F);
      }
      bBreakAmb <- F;

      nHomStart <- -1;
      nHomEnd <- -1;
      
    }
    nPrevState <- nState;
    nPrevMid <- nPosMid;
  }
}

