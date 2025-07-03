#setwd("/public3/group_crf/home/cuirf/bahaha_assembly/wgs/gatk/joined_genotype/snpeff_norm/rho")
source("asym.MK.R");
sOut <- "asymtotoic_mk_indels.txt"
arrBreaks <- c(seq(0.0, 0.8, 0.05) , 0.99, 1);
arrAlleleFreq <- arrBreaks[2:(length(arrBreaks)-1)]
AF_low <- 0.01
AF_high <- 0.99
dat <- read.table("consurf_af_sv.rec.bed", header=F, sep="\t")
dat <- dat[ , c(1:3, 6:ncol(dat) ) ]

dat$AF <- apply((dat[, 10:29]), 1, FUN= function(x) {nTotal = length(x)*2; p <- sum(x)/nTotal; return(p) })
dat$AFBin <- cut(dat$AF, breaks = arrBreaks);
dat <- dat[complete.cases(dat),];


#plot asymtotic MK-like :
dN <- sum(as.integer(dat$V10 == "HIGH" & dat$AF == 1));
dS <- sum(as.integer(dat$V10 == "LOW" & dat$AF == 1));
pN <- sum(as.integer(dat$V10 == "HIGH" & dat$AF < 1));
pS <- sum(as.integer(dat$V10 == "LOW" & dat$AF < 1));

arrPN <- c();
arrPS <- c();
arrAFBins <- as.character(levels(dat$AFBin));

for (sAFBin in arrAFBins[1:(length(arrAFBins)-1)]) {
  arrPN <- c(arrPN, sum(as.integer(dat$V10 == "HIGH" & dat$AFBin == sAFBin )));
  arrPS <- c(arrPS, sum(as.integer(dat$V10 == "LOW" & dat$AFBin == sAFBin )));
}

sink(file = sOut)
alpha <- 1 - ( (dS * pN) / (dN * pS) );
oConTable <- matrix(c(dN , dS, pN , pS ), nrow = 2);
oTest <- fisher.test(oConTable);

cat("dN\tdS\tpN\tpS\n");
cat(dN, "\t", dS, "\t", pN, "\t", pS, "\n");
cat("alpha global = ", alpha);
oTest

arrAsymMK <- asymptoticMK(d0=dS, d=dN, df=data.frame(f=arrAlleleFreq, p=arrPN, p0=arrPS, row.names=NULL), xlow=0.08, xhigh=0.9, output='table')
arrAsymMK

n10percIdx <- which(arrAlleleFreq == 0.1);
alpha_10perc <- 1 - ( (dS * arrPN[n10percIdx] ) / (dN * arrPS[n10percIdx] ) );
alpha_asymptotic <- arrAsymMK['alpha_asymptotic'];

alphaDiff <- alpha_asymptotic - alpha_10perc;

cat("alpha_10%\talpha_original\talpha_asymptotic\talphaDiff\n");
cat(alpha_10perc, "\t",  as.numeric(arrAsymMK['alpha_original']) , "\t", as.numeric(alpha_asymptotic), "\t", as.numeric(alphaDiff), "\n");


sink();
pdf(file=paste(sOut, ".pdf",sep="") , width=6, height=4)
asymptoticMK(d0=dS, d=dN, df=data.frame(f=arrAlleleFreq, p=arrPN, p0=arrPS, row.names=NULL), xlow=AF_low, xhigh=AF_high, output='default' )
abline(h=alpha_10perc , col='blue');
abline(h= as.numeric(arrAsymMK['alpha_original'])  , col='green');
abline(h= 0  , col='black');

dev.off();
  
