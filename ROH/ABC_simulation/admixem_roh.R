#setwd("D:\\data\\Ray\\working\\sysu\\server\\mop\\roh\\sim");

args = commandArgs(trailingOnly=TRUE)


sChrCoord <- args[1];
sGeno <- args[2];
nWinSize <- as.numeric(args[3]); # 50000;
nStepSize <- as.numeric(args[4]); # 10000;
nMaxInd <- as.numeric(args[5]); #20;#maximal sample 10 individuals

datCoord <- read.table(sChrCoord, header=T, stringsAsFactors = F)
datGeno <- read.table( sGeno , header=F, stringsAsFactors = F, nrows=nMaxInd*2);

nMarkers <- nrow(datCoord)
datGenotypes <- datGeno[, 7:(nMarkers+7-1)];


datHOM <- (datGenotypes[seq(1,nrow(datGenotypes),2) ,] == datGenotypes[seq(2,nrow(datGenotypes),2) ,])

arrWinSizes <- c();
nGenomeLen <- 0;
nInds <- nrow(datHOM);
if (nInds > nMaxInd) {
  nInds <- nMaxInd;
}
for(i in 1:nInds ) {
  datInd <- cbind(datCoord,datHOM[i,]);
  nLastMarkerCoord <- datCoord[nMarkers,1];
  nGenomeLen <- nGenomeLen + nLastMarkerCoord;
  datHetWin <- data.frame(start=numeric(), end=numeric(), mid=numeric(), het=numeric() ) ;
  nCounter <- 0;
  for (nWinStart in seq(1,nLastMarkerCoord-nWinSize/2, nStepSize) ) {
    nWinEnd <- nWinStart + nWinSize-1;
    nHet <- length( datInd[datInd[,1]>= nWinStart & datInd[,1]<= nWinEnd & (!datInd[,2]), 2]);
    nCounter <- nCounter +1;
    datHetWin[nCounter, ] <-  c(nWinStart, end=nWinEnd, mid=mean(nWinStart+nWinEnd)/2, het=nHet/nWinSize  );
  }
  
  nPrevState <- F;
  nPrevPos <- 0;
  nHomStart <- -1;
  for(nWin in 1:nrow(datHetWin)) {
    nHom <- (datHetWin[nWin, 'het']==0);
    nPos <- datHetWin[nWin, 'mid'];
    if (nHom && !nPrevState) {
      nHomStart <- nPos;
    }
    
    if ((!nHom) && nPrevState) {
      nWinLen <- nPrevPos - nHomStart + 1;
      arrWinSizes[length(arrWinSizes)+1] <- nWinLen;
      nHomStart <- -1
    }

    if (nHom && nPrevState && nWin==nrow(datHetWin)) {
      nWinLen <- nPrevPos - nHomStart + 1;
      arrWinSizes[length(arrWinSizes)+1] <- nWinLen;
      nHomStart <- -1
    }
    nPrevState <- nHom;
    nPrevPos <- nPos;
  }
  
}

arrWinSizes <- arrWinSizes[arrWinSizes>=nWinSize]
#hist(log10(arrWinSizes) );
if (length(arrWinSizes)>0) {
	nLogMean <- mean(log10(arrWinSizes));
	nLogSD <- sd(log10(arrWinSizes));
	nFROH <- sum(arrWinSizes)/nGenomeLen
	cat(nLogMean, "\t", nLogSD, "\t", nFROH)
} else {
	cat("0\t0\t0")
}
