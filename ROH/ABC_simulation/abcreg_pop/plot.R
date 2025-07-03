#setwd("/public3/group_crf/home/cuirf/macropodus_compare/roh_sims/abc/abcreg");
nQuant <- 0.025;
for (i in 0:2) {
sF <- paste0("reg.reseq.", i, ".tangent.post.gz");
oF <- gzfile(sF)

dat <- read.table(oF, header=F)
colnames(dat) <- c('gen', 'Ne', 'rec_per_meiosis');
hist(log10(dat$gen) )
hist(log10(dat$Ne) )
hist(log10(dat$rec_per_meiosis) )


hist(dat$gen )
hist(dat$Ne)
hist(dat$rec_per_meiosis )

datRet <- NULL;
nMeanGen <- mean(log10(dat$gen))
nSDGen <- sd(log10(dat$gen))
nLowerGen <- qnorm(nQuant, mean=nMeanGen, sd = nSDGen)
nUpperGen <- qnorm(1-nQuant, mean=nMeanGen, sd = nSDGen)
datRet <- rbind(datRet, data.frame(param="Log10_Gen", mean=nMeanGen, meanOrig=10^nMeanGen, sd=nSDGen, sdOrig=10^nSDGen,   lower=10^nLowerGen, upper=10^nUpperGen))

nMeanNe <- mean(log10(dat$Ne))
nSDNe <- sd(log10(dat$Ne))
nLowerNe <- qnorm(1e-4, mean=nMeanNe, sd = nSDNe)
nUpperNe <- qnorm(1-1e-4, mean=nMeanNe, sd = nSDNe)
datRet <- rbind(datRet, data.frame(param="Log10_Ne", mean=nMeanNe, meanOrig=10^nMeanNe, sd=nSDNe, sdOrig=10^nSDNe, lower=10^nLowerNe, upper=10^nUpperNe))


nMeanRec <- mean(log10(dat$rec_per_meiosis))
nSDRec <- sd(log10(dat$rec_per_meiosis))
nLowerRec <- qnorm(1e-4, mean=nMeanRec, sd = nSDRec)
nUpperRec <- qnorm(1-1e-4, mean=nMeanRec, sd = nSDRec)
datRet <- rbind(datRet, data.frame(param="Log10_Rec", mean=nMeanRec, meanOrig=10^nMeanRec, sd=nSDRec, sdOrig=10^nSDRec, lower=10^nLowerRec, upper=10^nUpperRec))

write.table(datRet, file=paste0(sF,".updated_prior.txt"), sep="\t", quote = F, row.names = F, col.names = T)
}
