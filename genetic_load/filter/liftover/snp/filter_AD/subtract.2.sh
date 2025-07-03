sF1=/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/liftover/allpops.cds.bothDP.bed
sF2=/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_bothSNP/snp/allpops.refMOP.NONSYN.discarded.bed

bedtools subtract -a $sF1 -b $sF2 > allpops.cds.bothDP.without_discarded.bed
