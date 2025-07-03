sF=/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_bothSNP/snp/allpops.refMHK.NONSYN.polarized.out.filtered.bed

bedtools intersect -a $sF -b allpops.merged.snpeff.vcf.AD.bed -wa -wb > allpops.refMHK.NONSYN.polarized.out.filtered.withAD.bed 
