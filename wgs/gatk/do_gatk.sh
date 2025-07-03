#!/bin/bash
#do gatk.sh for each sample.bam file in /data/projects/zwang/m.op/GATK/macro_for_add
for sSample in /data/projects/zwang/m.hk/GATK/final/*.bam; do
#	echo $sSample
        sBase=`basename $sSample`
        sName=${sBase/.bam/}
#	echo $sName
        cd /data/projects/zwang/m.hk/GATK/final/gatk/$sName
        source /data/projects/zwang/m.hk/GATK/final/gatk/gatk.sh $sSample &
done;
wait
