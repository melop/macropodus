module load bcftools
bcftools view -h joincalled_gvcf/joincalled.genotyped.g.vcf.gz | grep '^##contig=' | awk -F '[=,>]' '{print $3 "\t0\t" $5"\t0\t" $5"\t0\t" $5}' | head -23 > alltabs.tab
