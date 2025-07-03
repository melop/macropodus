#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/slurm-%a.out
#SBATCH -e slurmlog/slurm-%a.err
#SBATCH -a 0-43
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8


module load clusterbasics
module load rohan
ref="/public3/group_crf/home/g20wangzhx36/m.hk/nextdenovo/mhk.LG.nextdenovo.sgs.2.sorted.fa"
sBam_path="/fast3/group_crf/home/g20wangzhx36/m.hk/GATK/final2/merge_samples/ROHan/all_dedup_bam"
nJobCount=0
for i in $sBam_path/*.filtered.bam; do
	if (( ${SLURM_ARRAY_TASK_ID} == $nJobCount  )); then
	sF=$(basename $i)
	sStem=${sF/.filtered.bam/}
	echo $sStem
	( mkdir -p $sStem && cd $sStem ; \
	rohan --rohmu 2e-5 -t 8 --size 50000 --step 100 -o $sStem $ref $i > run.log 2>&1 ; cd .. )
	fi
	nJobCount=$(( nJobCount + 1 ))
done;


