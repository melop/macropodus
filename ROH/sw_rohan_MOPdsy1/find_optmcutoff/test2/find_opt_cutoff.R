setwd("/data2/projects/zwang/m.hk/ROH/simulation1/find_optmcutoff");

arrCutoffs <- seq(1e-5, 0.1, 1e-5)

datROHAN <- read.table("X5.07.rohan.het.bed", sep="\t", header=F, stringsAsFactors = F);
datROHAN <- datROHAN[complete.cases(datROHAN),];
datROHAN <- datROHAN[datROHAN$V1 %in% paste0("mhkscf_", 1:23) , ];

colnames(datROHAN) <- c('chr', 'winstart', 'winend', 'esthet');
datROHAN$coord <- paste(datROHAN$chr, datROHAN$winstart, datROHAN$winend, sep=":")

datROHWin <- read.table("X5.07.roh.wins.bed" , sep="\t", header=F);
datROHWin$coord <- paste(datROHWin$V1, datROHWin$V2, datROHWin$V3, sep=":")

datROHWin$ROH <- T;

datROHAN$ROH <-F
datROHAN$ROH[datROHAN$coord %in% datROHWin$coord ] <- T;

nFROH <- sum(datROHAN$ROH)/nrow(datROHAN)


datSimRets <- NULL

for(nCutoff in arrCutoffs) {
  nROH_on_ROH <- sum((datROHAN$esthet < nCutoff) & (datROHAN$ROH));
  nROH_on_HET <- sum((datROHAN$esthet < nCutoff) & (datROHAN$ROH==F));
  nHET_on_ROH <- sum((datROHAN$esthet >= nCutoff) & (datROHAN$ROH));
  nHET_on_HET <- sum((datROHAN$esthet >= nCutoff) & (datROHAN$ROH==F));
  nEstFROH <- (nROH_on_ROH + nROH_on_HET) / (nROH_on_ROH + nROH_on_HET + nHET_on_ROH + nHET_on_HET);
  datSimRets <- rbind(datSimRets, data.frame(Cutoff = nCutoff,nROH_on_ROH=nROH_on_ROH, nROH_on_HET=nROH_on_HET, nHET_on_ROH=nHET_on_ROH, nHET_on_HET=nHET_on_HET, estroh = nEstFROH ) );
}

datSimRets$correct <- (datSimRets$nROH_on_ROH + datSimRets$nHET_on_HET) / rowSums(datSimRets[,2:5])
plot(datSimRets$Cutoff , datSimRets$correct, type="l")

# 找到最大y值对应的索引
max_correct_idx <- which.max(datSimRets$correct)
# 获取对应的x值（Cutoff）
x_value_at_max_y <- datSimRets$Cutoff[max_correct_idx]
print(x_value_at_max_y)
