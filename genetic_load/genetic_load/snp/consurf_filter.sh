sF=/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_bothSNP/snp/allpops.refMHK.NONSYN.polarized.out.filtered.bed

bedtools intersect -a $sF -b consurf_scores_flat.bed -loj > allpops.refMHK.NONSYN.polarized.out.filtered.withconsurf.bed
