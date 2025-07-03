MOP=/data2/projects/zwang/m.op/Snpeff/all_MHK_MOP/allpops.withoutgroup.g.vcf.gz
bcftools query -f '%CHROM\t%POS0\t%END\t[%AD\t]\n' $MOP > test.allpops.refMOP.gvcf.AD.bed
