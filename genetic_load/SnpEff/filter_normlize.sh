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
genome=/public3/group_crf/home/g20wangzhx36/m.hk/release1.0/genome/scf.softmasked.fa

#bcftools view -O z -M 2 -m 2 -i 'TYPE=="snp" || TYPE=="indel"' allpops.withoutgroup.g.vcf.gz > allpops.withoutgroup.vcf.gz #run in albonubes

nJobCount=0
for i in $(bcftools view -h allpops.withoutgroup.vcf.gz | grep "^##contig=<ID=" | cut -f3 -d'=' | cut -f1 -d','); do
	if (( ${SLURM_ARRAY_TASK_ID} == $nJobCount  )); then
#		echo $i
		( bcftools norm -w 99999999 -O z -f $genome -r $i allpops.withoutgroup.vcf.gz > allpops.withoutgroup.norm.$i.vcf.gz
		tabix allpops.withoutgroup.norm.$i.vcf.gz )
	fi
	nJobCount=$(( nJobCount + 1 ))
done;

#bcftools merge -o allpops.withoutgroup.norm.vcf.gz -O z allpops.withoutgroup.norm.*.vcf.gz
#tabix allpops.withoutgroup.norm.vcf.gz

