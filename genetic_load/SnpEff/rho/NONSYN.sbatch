#!/usr/bin/bash
#SBATCH -p gpu,himem,hugemem,blade
#SBATCH -o slurmlog/slurm-%a.out
#SBATCH -e slurmlog/slurm-%a.err
#SBATCH -a 0-542
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load bcftools
module load samtools

REF=Macropodus_hongkongensis

nJobCount=0
for i in $(bcftools view -h allpops.withoutgroup.vcf.gz | grep "^##contig=<ID=" | cut -f3 -d'=' | cut -f1 -d','); do
	if (( ${SLURM_ARRAY_TASK_ID} == $nJobCount  )); then
#		echo $i
		( bcftools view -r $i allpops.merged.snpeff.renameoutgroup.vcf.gz -o allpops.snpeff.renameoutgroup.$i.vcf.gz
		tabix allpops.snpeff.renameoutgroup.$i.vcf.gz
		SnpEffVcf=allpops.snpeff.renameoutgroup.$i.vcf.gz
		Out=NONSYN.$i.polarized.out.txt
       		php getNonsyn.php "$SnpEffVcf" "$Out")
	fi
	nJobCount=$(( nJobCount + 1 ))
done;


