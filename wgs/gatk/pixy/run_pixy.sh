#!/usr/bin/bash
#SBATCH -p blade,gpu
#SBATCH -c 32

module load bcftools
module load samtools
module load vcftools
module load R
source /public/apps/miniconda3/etc/profile.d/conda.sh

#less DP.txt | awk 'BEGIN{n=0;min=0};{if ($1!="FCG:" && $1!="YN-mop6:" && $2=="lower"){n+=1;min+=$4}};END{print min/n}'


##filter
vcftools --gzvcf all.samples_refmhk.g.vcf.gz \
	--remove-indels \
	--max-missing 0.8 \
	--min-meanDP 4 \
	--max-meanDP 44 \
	--recode --stdout | bgzip > all.samples_refmhk.filtered.vcf.gz;
#gzip -d all.samples_refmhk.filtered.vcf.gz;
#bgzip -f all_rmFCG_yn6.filtered.g.vcf;
tabix all.samples_refmhk.filtered.vcf.gz;

#conda activate pixy and run pixy
conda activate pixy
pixy --stats pi fst dxy \
	--vcf all.samples_refmhk.filtered.vcf.gz \
	--populations samples.txt \
	--window_size 10000 \
	--n_cores 32 \
	--output_folder output \
	--output_prefix pixy_output
