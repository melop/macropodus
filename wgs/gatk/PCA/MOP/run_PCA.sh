plink=/data/software/Plink/plink
Species=MOP

# Generate the input file in plink format
$plink --vcf all.mop_selected_refmhk_filtersnps.vcf.gz --make-bed --out $Species --allow-extra-chr --const-fid ;


$plink --allow-extra-chr --threads 20 -bfile $Species --pca 20 --out $Species;
