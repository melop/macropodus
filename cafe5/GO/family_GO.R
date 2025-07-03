setwd("/data2/projects/zwang/macropodus_compare/cafe5_RC/GO/with_genesymbol");
library(org.Dr.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)


sInDir <- "../../runresults_k4/7/";
sRefSp <- sSp <- "MacropodusOpercularis";
#sRefSp <- sSp <- "MacropodusHongkongensis"
nFDRCutoff <- 0.2
sCategory <- 'BP' ; #biological process
datPGID2Ref <- read.table("/data/projects/zwang/macropodus_compare/synteny/genespace/4spp/rundir/orthofinder/Results_Nov30/Orthogroups/Orthogroups.tsv", header = T, stringsAsFactors = F, fill=T, sep="\t")
datFam2Gene <- read.table("/data2/projects/zwang/macropodus_compare/synteny/genespace_RC/rundir/orthofinder/Results_Nov01/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", header=T, sep="\t", stringsAsFactors = F);
datFam2Gene$OG <- gsub('N0.', '', datFam2Gene$OG , fixed = T)

datChange <- read.table(paste0(sInDir, "/Gamma_change.tab" ), sep="\t", header=T);
datChange <- datChange[, c(1, grep(sSp, colnames(datChange), fixed = T)) ];

datBranchProb <- read.table(paste0(sInDir, "/Gamma_branch_probabilities.tab" ), sep="\t", header=T, comment.char = '/', fill = T);

datChange <- merge(datChange , datBranchProb[, c(1,grep(sSp, colnames(datBranchProb)))], by.x = "FamilyID", by.y=1)

datChange$fdr <- p.adjust(datChange[, 3])

arrExpandedFam <- datChange[datChange$fdr<=nFDRCutoff & datChange[,2] > 0  ,1]
arrShrinkedFam <- datChange[datChange$fdr<=nFDRCutoff & datChange[,2] < 0  ,1]
arrControlFam <- datChange[ ,1]

fnFam2Gene <- function(arr) {
  arrGenes <- datFam2Gene[datFam2Gene$OG %in% arr, sSp];
  return(unique(str_trim(unlist(strsplit(x = arrGenes, split = ',')))));
}

arrExpandedGenes <- fnFam2Gene(arrExpandedFam);
arrShrinkedGenes <- fnFam2Gene(arrShrinkedFam);
arrControlGenes <- fnFam2Gene(arrControlFam);

#prepare mapping tables

datRNAID2ZebraFishHuman <- datPGID2Ref[, c(sRefSp,  'Human',  'Zebrafish') ];
datRNAID2ZebraFishHuman <- datRNAID2ZebraFishHuman[(!is.na(datRNAID2ZebraFishHuman[,sRefSp])) & (!datRNAID2ZebraFishHuman[,sRefSp]=='') , ];
datRNAID2ZebraFishHuman <- datRNAID2ZebraFishHuman[!(datRNAID2ZebraFishHuman[, 'Human']=='' & datRNAID2ZebraFishHuman[, 'Zebrafish']=='' ), ];


datZebrafishTrans2GeneID <- read.table("/data/projects/rcui/bahaha_assembly/wgs/sex_chrm_id_q0q1/GO/zebrafish_ncbi_ensembl_id_map.txt", sep="\t", header=F) ;#as.data.frame( org.Dr.egENSEMBLTRANS)
datZebrafishTrans2GeneID$V5 <- gsub('\\.[0-9]*' , "", datZebrafishTrans2GeneID$V5)
datZebrafishTrans2GeneID <- datZebrafishTrans2GeneID[, c(2,5)]
colnames(datZebrafishTrans2GeneID) <- c("ZebrafishGeneID", 'trans_id');

datHumanTrans2GeneID <- read.table("/data/projects/rcui/bahaha_assembly/wgs/sex_chrm_id_q0q1/GO/human_ncbi_ensembl_id_map.txt", sep="\t", header=F); #as.data.frame( org.Hs.egENSEMBLTRANS)
datHumanTrans2GeneID$V5 <- gsub('\\.[0-9]*' , "", datHumanTrans2GeneID$V5)
datHumanTrans2GeneID <- datHumanTrans2GeneID[, c(2,5)]
colnames(datHumanTrans2GeneID) <- c("HumanGeneID", 'trans_id');


datHumanGOMap <- as.data.frame(org.Hs.egGO)
datZebrafishGOMap <- as.data.frame(org.Dr.egGO)

datMap <- NULL;
sTERM2GO <- paste0(sSp,"_term2go.tsv");
bAppend <- F;

sink(file=paste0(sSp, ".gomap.log"))
for(i in 1:nrow(datRNAID2ZebraFishHuman)) {
  arrZebrafishTransID <- str_trim(unlist(strsplit(x = datRNAID2ZebraFishHuman[i, 'Zebrafish'], split = ',')));
  arrHumanTransID <- str_trim(unlist(strsplit(x = datRNAID2ZebraFishHuman[i, 'Human'], split = ',')));
  arrSpTransID <- str_trim(unlist(strsplit(x = datRNAID2ZebraFishHuman[i, 1], split = ',')));
  arrZebrafishGeneID <- datZebrafishTrans2GeneID[datZebrafishTrans2GeneID$trans_id %in% arrZebrafishTransID, 'ZebrafishGeneID']
  arrHumanGeneID <- datHumanTrans2GeneID[datHumanTrans2GeneID$trans_id %in% arrHumanTransID, 'HumanGeneID']
  
  arrGOs <- datZebrafishGOMap[datZebrafishGOMap$gene_id %in% arrZebrafishGeneID & datZebrafishGOMap$Ontology ==sCategory, 'go_id']
  arrGOs <- c(arrGOs, datHumanGOMap[datHumanGOMap$gene_id %in% arrHumanGeneID & datHumanGOMap$Ontology ==sCategory, 'go_id']);
  arrGOs <- unique(arrGOs)
  if (length(arrGOs) == 0) {
    cat("GO not found: ", datRNAID2ZebraFishHuman[i, 1], "\n");
    next;
  } else {
    cat("GO found: ", datRNAID2ZebraFishHuman[i, 1], "\n");
  }
  for(sSpTransID in arrSpTransID) {
    datToBind <- data.frame(term=arrGOs, gene=sSpTransID, stringsAsFactors = F);
    write.table(datToBind, file = sTERM2GO, append = bAppend, sep = "\t", col.names = (!bAppend), row.names = F, quote = F)
    bAppend<-T;
    
  }
}
sink()

datMap <- read.table(sTERM2GO, sep="\t", header=T, stringsAsFactors = T);
#-------------------------------------------------------------------------------scripts above get GO term mapping to MOP gene
datRNAID2ZebraFishHuman$Human <- sub(",.*", "", datRNAID2ZebraFishHuman$Human)
datRNAID2ZebraFishHuman$Zebrafish <- sub(",.*", "", datRNAID2ZebraFishHuman$Zebrafish)

colnames(datRNAID2ZebraFishHuman)[1] <- "Gene";
datRNAID2ZebraFishHuman$Gene <- gsub('.', '_', datRNAID2ZebraFishHuman$Gene, fixed = T);

datOnlyCompHighImpact <- datRNAID2ZebraFishHuman;
datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datZebrafishTrans2GeneID, by.x="Zebrafish", by.y = "trans_id", all.x =T, all.y =F);
datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datHumanTrans2GeneID, by.x="Human", by.y = "trans_id", all.x =T, all.y =F);

#get gene symbol map
datHumanSymbolMap <- as.data.frame(org.Hs.egSYMBOL)
datZebrafishSymbolMap <- as.data.frame(org.Dr.egSYMBOL)

datSymbol <- datOnlyCompHighImpact[, c('Gene', "ZebrafishGeneID")]
datSymbol <- datSymbol[complete.cases(datSymbol) , ];
datSymbol <- merge(datSymbol, datZebrafishSymbolMap[,1:2], by.y="gene_id", by.x="ZebrafishGeneID", all.x=T, all.y=F )
datSymbol <- datSymbol[complete.cases(datSymbol) , c(3,2)];

datSymbolHs <- datOnlyCompHighImpact[, c('Gene', "HumanGeneID")]
datSymbolHs <- datSymbolHs[complete.cases(datSymbolHs) , ];
datSymbolHs <- merge(datSymbolHs, datHumanSymbolMap[,1:2], by.y="gene_id", by.x="HumanGeneID", all.x=T, all.y=F )
datSymbolHs <- datSymbolHs[complete.cases(datSymbolHs) , c(3,2)];

datSymbol <- rbind(datSymbol , datSymbolHs[ !(datSymbolHs$Gene %in% datSymbol$Gene),] );
datSymbol <- datSymbol[!duplicated(datSymbol),];
#-------------------------------------------------------------------------------script above get gene symbol mapping to MOP gene
#fnAddGeneSymbols <- function(datOut) {
#  arrSymbolStr <- c();
#  if (nrow(datOut) ==0) {
#    return(datOut)
#  }
#  for(nRow in 1:nrow(datOut)) {
#    
#    arrPGIDs <- unlist(strsplit(datOut$geneID[nRow] ,'/'));
#    arrGeneSymbols <- datSymbol[datSymbol$Gene %in% arrPGIDs, 1];
#    sStr <- "";
#    if (length(arrGeneSymbols)>0) {
#      sStr <- paste(arrGeneSymbols, collapse = '/');
#    }
#    arrSymbolStr <- c(arrSymbolStr, sStr);
#  }
#  datOut$genesymbols <- arrSymbolStr;
#  return(datOut)
#}
fnAddGeneSymbols <- function(datOut, datSymbol) {
  arrSymbolStr <- c()
  if (nrow(datOut) == 0) {
    return(datOut)
  }
  for (nRow in 1:nrow(datOut)) {
    arrPGIDs <- unlist(strsplit(datOut$geneID[nRow], '/'))
    sSymbolStr <- ""
    for (geneSymbolRow in 1:nrow(datSymbol)) {
      # 先拆分datSymbol中Gene列的每个元素（以'/'分割）
      arrGenesInSymbol <- unlist(strsplit(datSymbol[geneSymbolRow, "Gene"], ', ')) #逗号+空格
      matchFlag <- FALSE
      for (geneInSymbol in arrGenesInSymbol) {
        if (geneInSymbol %in% arrPGIDs) {
          matchFlag <- TRUE
          break
        }
      }
      if (matchFlag) {
        symbolValue <- datSymbol[geneSymbolRow, 1]
        sSymbolStr <- if (sSymbolStr == "") {
          symbolValue
        } else {
          paste(sSymbolStr, symbolValue, sep = "/")
        }
      }
    }
    arrSymbolStr <- c(arrSymbolStr, sSymbolStr)
  }
  datOut$genesymbols <- arrSymbolStr
  return(datOut)
}

#enrichment test:

# arrExpandedGenes 
# arrShrinkedGenes 
# arrControlGenes

sum(as.integer(unique(datMap$gene) %in% arrExpandedGenes))
sum(as.integer(unique(datMap$gene) %in% arrControlGenes))

#View( go2term(unique(datMap$term[datMap$gene %in% arrExpandedGenes])) )
arrExpandedGenes[!(arrExpandedGenes %in% unique(datMap$gene))]
                                                
oEnrichRet <- enricher(
  gene=arrExpandedGenes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = arrControlGenes ,
  minGSSize = 10,
  maxGSSize = 5000,
  qvalueCutoff = 0.1,
  TERM2GENE = datMap ,
  TERM2NAME = go2term(unique(datMap$term) )
)
datOut <- oEnrichRet@result[ oEnrichRet@result$p.adjust <=1, ]
datOut$geneID <- gsub('.', '_', datOut$geneID, fixed = T);
#View(datOut)
write.table(fnAddGeneSymbols(datOut,datSymbol), paste0(sSp, "_expandedfam.GO.with_genesymbol.tsv"), sep = "\t", col.names = T, row.names = F, quote = F);


oEnrichRet <- enricher(
  gene=arrShrinkedGenes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = arrControlGenes ,
  minGSSize = 10,
  maxGSSize = 5000,
  qvalueCutoff = 0.1,
  TERM2GENE = datMap ,
  TERM2NAME = go2term(unique(datMap$term) )
)
datOut <- oEnrichRet@result[ oEnrichRet@result$p.adjust <=1, ]
datOut$geneID <- gsub('.', '_', datOut$geneID, fixed = T);
#View(datOut)
write.table(fnAddGeneSymbols(datOut,datSymbol), paste0(sSp, "_shrinkedfam.GO.with_genesymbol.tsv"), sep = "\t", col.names = T, row.names = F, quote = F);

arrUnMappedExpanded <- arrExpandedGenes[!(arrExpandedGenes %in% unique(datMap$gene))];

arrShow <-c();
for(sUnmapped in arrUnMappedExpanded) {
  arrShow <- c(arrShow, grep(sUnmapped, datPGID2Ref$BahahaTaipingensis));
}

#View(datPGID2Ref[unique(arrShow),]);
