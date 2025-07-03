#get a table of variant stats from chromsome1 for visualization 
#in order to pick the appropriate cutoffs for variant filtering
chrom=mhkscf_1

echo -e 'POS\tGT\tAD\tDP\tQD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum\n' > var.stats.$chrom.tsv
bcftools query --include 'type="snp" && CHROM="'$chrom'"' --format '[%POS]\t[%GT]\t[%AD]\t%DP\t%QD\t%FS\t%SOR\t%MQ\t%MQRankSum\t%ReadPosRankSum\n' DSY-1.genotyped.g.vcf.gz >> var.stats.$chrom.tsv
