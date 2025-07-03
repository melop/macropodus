plink=/data/software/Plink/plink
Name=all_refmhk_new
#STEP 1
bcftools view -O z -o all.samples_refmhk.filter_forplink.snps.vcf.gz -i 'count(REF)==1 && count(ALT)==1 && QUAL>=60 && TYPE=="snp" && INFO/DP > 462 && INFO/DP < 3716 && FILTER!="LowQual" && F_PASS(GT!="mis")>0.90' all.samples_refmhk.g.vcf.gz

#STEP2
# Generate the input file in plink format
$plink --vcf all.samples_refmhk.filter_forplink.snps.vcf.gz --make-bed --out $Name --allow-extra-chr --const-fid 0

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' $Name.bim > $Name.bim.tmp 
mv $Name.bim.tmp $Name.bim

bash run_admixture.sh > run_admixture.log 2>&1 
