MHK=/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/allpops.withoutgroup.g.vcf.gz
bcftools query -f '%CHROM\t%POS0\t%END\t[%AD\t]\n' $MHK > test.allpops.refMHK.gvcf.AD.bed
