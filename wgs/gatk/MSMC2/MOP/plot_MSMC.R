#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/PLPv1.2LG/msmc2");
setwd("/data2/projects/zwang/m.hk/MSMC_suppl/MOP")
pdf("MOP.MSMC.noBS.pdf", width=12, height=8);
#pdf("praslin.pdf", width=7, height=6);
#pdf("mahesouth.pdf", width=7, height=6);
#pdf("allwild.pdf", width=7, height=6);

sMainFolder <- "msmc2ret/";
sBSFolder <- "bootstrapped/"
nBSReps <- 30;
arrYLim <- c(1e2, 1e9);
arrXLim <- c(1e3, 1e7);

nMu <- 2.7E-9; # estimated from divergence between BTP and Miichthys, this is per generation. assuming 10 years of generation time, per year is 7.75273E-10


fnPlotPop <- function(sPop, arrCol, bBSPlot = F, bAppendPlot=T, sThisBSFile1 = "") {
 
  sInFile1 <- paste(sMainFolder,"/",sPop, ".final.txt", sep="");
  arrCol <- add.alpha(arrCol, 0.8);
  
  if (bBSPlot) {
    cat("bs: ", sThisBSFile1,"\n");
    sInFile1 <- sThisBSFile1;
    arrCol <- add.alpha(arrCol, 0.1);
  }
  
  cat("Open ", sInFile1,"\n");
  datMSMC1 <- read.table(sInFile1, header=T, sep="\t" );
  datMSMC1$gen <- datMSMC1$left_time_boundary/nMu;
  datMSMC1$gen[1] <- 0.01;
  datMSMC1$popsize <- (1/datMSMC1$lambda)/(2*nMu);
  
  
  if (bBSPlot == F && (!bAppendPlot) ) {
    plot( datMSMC1$gen ,datMSMC1$popsize ,  type = 's' , log='xy', xlim=arrXLim, ylim = arrYLim, xlab="Generations" , ylab = "Pop size", lwd=3, col=arrCol[1])
    # 设置当前图形的坐标轴字体大小
    axis(1, cex.axis = 1.8)
    axis(2, cex.axis = 1.8)
    mtext("Generations", side = 1, line = 3, cex = 2.4)
    mtext("Pop size", side = 2, line = 3, cex = 2.4)
  } else {
    lines( datMSMC1$gen ,datMSMC1$popsize ,  type = 's' , lwd=3, col=arrCol[1])
    
  }
  

  if (bBSPlot==TRUE) {
    #plot variation
    return();
  } else {
    cat("Plot bs...\n");
    for (nRep in 1:nBSReps) {
      sBSF1 <- paste(sBSFolder,'/',sPop, '/_', nRep, '/out.final.txt' , sep="");
      #cat(sBSF1);
      if (file.exists(sBSF1)) {
        fnPlotPop(sPop, arrCol, bBSPlot = T, sThisBSFile1=sBSF1 )
      }
    }
    
    # datReal <- read.table( paste(sRealFolder, '/',sRealFile ,sep=""), header=T, sep="\t");
    # points(datReal$Pop1_Gen, datReal$Pop1_Popsize, col=arrCol[2])
    # points(datReal$Pop2_Gen, datReal$Pop2_Popsize, col=arrCol[3])
    # 
  }
  
  #construct output table:
  return();
  
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# 创建 7 个容易区分的暖色调颜色
colors <- c("#2873B3","#2EBEBE","#74B346", "#F1CC2F","#7D4444","#A14462","#8264CC")

#Wild Praslin:

datSamples <- read.table("samples.txt", sep=' ', header=F);
arrSamples <- datSamples$V1;
#arrCol <- rainbow(length(arrSamples));
bAppend <-F;
# 创建一个空的向量来存储种群名称和颜色的对应关系
legend_text <- character()

for(i in 1:length(arrSamples) ) {
	fnPlotPop(arrSamples[i], colors[i] , bAppendPlot = bAppend);
	bAppend <-T
	# 添加注释文本
	legend_text <- c(legend_text, arrSamples[i])
}

# 添加图例
legend("topright", legend = legend_text, col = colors, lwd = 3, bty = "n", cex = 1.8)

dev.off()
#----------------------------------------------------------------------------------------



