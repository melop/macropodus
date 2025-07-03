CPU=3
REF=/data/projects/zwang/m.hk/GATK/mhk.LG.nextdenovo.sgs.2.sorted.fa
for sR1 in /data/projects/zwang/m.hk/GATK/final/trimmed_macro/*.paired_1.fq.gz; do
	sDir=`dirname $sR1`
	sBase=`basename $sR1`
	sSample=${sBase/_raw.paired_1.fq.gz/}
	sR2=${sDir}/${sSample}_raw.paired_2.fq.gz
	( bwa mem -t $CPU $REF $sR1 $sR2 \
        | samtools view  -u - | samtools sort - -m 10g -o $sSample.bam ) > $sSample.log 2>&1 &
#       echo $sR1
#       echo $sDir
#	echo $sBase
#	echo $sSample
#	echo $sR2
done
wait
