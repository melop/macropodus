allpops=/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_DP/refMHK/allpops.g.vcf.gz
sPath=/data2/projects/zwang/m.hk/all_gvcf/
bcftools merge -o allpops.withoutgroup.g.vcf.gz -O z --threads 32 $allpops $sPath/YN-12*.gz $sPath/WL-1*.gz
tabix allpops.withoutgroup.g.vcf.gz
