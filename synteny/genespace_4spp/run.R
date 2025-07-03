#setwd("~/bahaha_assembly/synteny/genespace/example/");
library(GENESPACE)
runwd <- file.path("/data/projects/zwang/macropodus_compare/synteny/genespace/4spp/rundir")
#list.files(runwd, recursive = T, full.names = F)

gpar <- init_genespace(
  genomeIDs = c('MacropodusOpercularis', 'MacropodusHongkongensis', 'Zebrafish', 'Human'),
  speciesIDs = c('Macropodus_opercularis', 'Macropodus_hongkongensis', 'Zebrafish', 'Human' ),
  versionIDs = c("1.0", "1.0", "107", "107"),
  ploidy = rep(1,4),
  diamondMode = "default",
  orthofinderMethod = "default",
  wd = runwd,
  nCores = 80,
  minPepLen = 30,
  gffString = "gff",
  pepString = "pep",
  path2orthofinder = "/opt/miniconda3/bin/orthofinder",
  path2diamond = "/opt/miniconda3/bin/diamond",
  path2mcscanx = "/data/software/MCScanX/",
  rawGenomeDir = file.path(runwd, "rawGenomes"))


parse_annotations(
  gsParam = gpar,
  gffEntryType = "gene",
  gffIdColumn = "locus",
  gffStripText = "locus=",
  headerEntryIndex = 1,
  headerSep = " ",
  headerStripText = "locus=")

gpar <- run_orthofinder(
  gsParam = gpar)

gpar <- synteny(gsParam = gpar)

arrCol <- rainbow(23);
regs <- data.frame(
  genome = rep('MacropodusOpercularis', 23),
  chr =1:23)

datInvert <- data.frame(genome="MacropodusHongkongensis", chr=c(1, 4, 8, 10, 6, 18,19,20, 22) )

ripdat <- plot_riparianHits(
  gpar,onlyTheseRegions = regs, refChrCols = arrCol,minGenes2plot=50, invertTheseChrs = datInvert,
  blackBg = F, chrFill = "orange",returnSourceData=T,
  chrBorder = "grey", useOrder=F, labelTheseGenomes = c('MacropodusOpercularis', 'MacropodusHongkongensis') )

#dump('ripdat', file="ripdata.R");

pg <- pangenome(gpar)
