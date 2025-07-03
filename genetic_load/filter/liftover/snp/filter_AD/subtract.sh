MHK=/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/allpops.withoutgroup.g.vcf.gz
MOP=/data2/projects/zwang/m.op/Snpeff/all_MHK_MOP/allpops.withoutgroup.g.vcf.gz
bed=/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/liftover/snp/allpops.cds.bothDP.without_discarded.bed

#awk '{print $1"\t"$2"\t"$3}' $bed > temp.refMOP.bed
#awk '{print $4"\t"$5"\t"$6}' $bed > temp.refMHK.bed

bcftools view -R temp.refMOP.bed $MOP -Oz -o allpops.cds.bothDP.without_discarded.refMOP.vcf.gz &

bcftools view -R temp.refMHK.bed $MHK -Oz -o allpops.cds.bothDP.without_discarded.refMHK.vcf.gz &
wait
