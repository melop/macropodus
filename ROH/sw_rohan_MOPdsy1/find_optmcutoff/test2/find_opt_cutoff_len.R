#setwd("/data/projects/rcui/mhk/nextdenovo/gatk_vcf/slidingwin_het/set_cutoff2");
setwd("/data2/projects/zwang/m.hk/ROH/sw_rohan_DSY1/find_optmcutoff/test2")
arrCutoffs <- seq(1e-5, 0.1, 1e-5)
datROHTruth <- read.table("roh.bed", sep="\t", header=F);
datROHTruth$len <- datROHTruth$V3 - datROHTruth$V2

datROHAN <- read.table("X22.rohan.het.bed", sep="\t", header=F, stringsAsFactors = F);
datROHAN <- datROHAN[complete.cases(datROHAN),];
datROHAN <- datROHAN[datROHAN$V1 %in% unique(datROHTruth$V1) , ];

colnames(datROHAN) <- c('chr', 'winstart', 'winend', 'esthet');

datROHWin <- read.table("X22.roh.wins.bed" , sep="\t", header=F);
datROHWin <- datROHWin[complete.cases(datROHWin),]

arrChrs <- unique(datROHAN$chr);

datSimRets <- NULL

for(nCutoff in arrCutoffs) {
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

 # }


  nEstFROH <- (nROH_on_ROH + nROH_on_HET) / (nROH_on_ROH + nROH_on_HET + nHET_on_ROH + nHET_on_HET);
  datSimRets <- rbind(datSimRets, data.frame(Cutoff = nCutoff,nROH_on_ROH=nROH_on_ROH, nROH_on_HET=nROH_on_HET, nHET_on_ROH=nHET_on_ROH, nHET_on_HET=nHET_on_HET, estroh = nEstFROH ) );
}

datSimRets$correct <- (datSimRets$nROH_on_ROH + datSimRets$nHET_on_HET) / rowSums(datSimRets[,2:5])

datSimRets$TP <- datSimRets$nROH_on_ROH  / rowSums(datSimRets[,c(2,4)])
datSimRets$TN <- datSimRets$nHET_on_HET / rowSums(datSimRets[,c(3,5)])
datSimRets$FP <- datSimRets$nROH_on_HET/ rowSums(datSimRets[,c(3,5)])
datSimRets$FN <- datSimRets$nHET_on_ROH/ rowSums(datSimRets[,c(2,4)])

plot(datSimRets$Cutoff , datSimRets$correct, type="l", main="correct", xlim=c(0,0.01));
abline(v=0.00082, col="red")

plot(datSimRets$Cutoff , datSimRets$TP, type="l", main="TP", xlim=c(0,0.01));
abline(v=0.00082, col="red")
plot(datSimRets$Cutoff , datSimRets$TN, type="l", main="TN", xlim=c(0,0.01));
abline(v=0.00082, col="red")
plot(datSimRets$Cutoff , datSimRets$FP, type="l", main="FP", xlim=c(0,0.01));
abline(v=0.00082, col="red")
plot(datSimRets$Cutoff , datSimRets$FN, type="l", main="FN", xlim=c(0,0.01));
abline(v=0.00082, col="red")

