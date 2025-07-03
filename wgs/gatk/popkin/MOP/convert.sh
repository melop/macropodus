bcftools view  -M2 -m2 -v snps all.mop_selected_refmhk_filtersnps.vcf.gz -O z > in.vcf.gz 
/data/software/plink2 --allow-extra-chr --vcf in.vcf.gz  --make-bed --out converted
