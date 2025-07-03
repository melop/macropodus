setwd("/public3/group_crf/home/cuirf/macropodus_compare/roh_sims/abc/abcreg");

sF <- "reg.reseq.0.tangent.post.gz"
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
nLowerGen <- qnorm(1e-4, mean=nMeanGen, sd = nSDGen)
nUpperGen <- qnorm(1-1e-4, mean=nMeanGen, sd = nSDGen)
datRet <- rbind(datRet, data.frame(param="Log10_Gen", mean=nMeanGen, sd=nSDGen, lower=nLowerGen, upper=nUpperGen))

nMeanNe <- mean(log10(dat$Ne))
nSDNe <- sd(log10(dat$Ne))
nLowerNe <- qnorm(1e-4, mean=nMeanNe, sd = nSDNe)
nUpperNe <- qnorm(1-1e-4, mean=nMeanNe, sd = nSDNe)
datRet <- rbind(datRet, data.frame(param="Log10_Ne", mean=nMeanNe, sd=nSDNe, lower=nLowerNe, upper=nUpperNe))


nMeanRec <- mean(log10(dat$rec_per_meiosis))
nSDRec <- sd(log10(dat$rec_per_meiosis))
nLowerRec <- qnorm(1e-4, mean=nMeanRec, sd = nSDRec)
nUpperRec <- qnorm(1-1e-4, mean=nMeanRec, sd = nSDRec)
datRet <- rbind(datRet, data.frame(param="Log10_Rec", mean=nMeanRec, sd=nSDRec, lower=nLowerRec, upper=nUpperRec))

write.table(datRet, file=paste0(sF,".updated_prior.txt"), sep="\t", quote = F, row.names = F, col.names = T)
