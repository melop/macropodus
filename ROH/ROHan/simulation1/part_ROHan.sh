picard="/data/software/picard.jar"
rohan="/data/software/ROHan-1.0.1/bin/rohan"
ref="/data/projects/rcui/mhk/nextdenovo/nextpolish_stlfr2/mhk.LG.nextdenovo.sgs.2.sorted.fa"

sIn="/data/projects/zwang/m.hk/GATK/final/merge_samples/ROHan/all_dedup_bam/SZ-1.filtered.bam"
sBam_path="/data2/projects/zwang/m.hk/ROH/simulation1.2"
arr=(2 4 6 8 10 12 14 16 18 20 22 24)
for i in ${arr[@]}; do
#	a=`echo "scale=4;$i/185.07"|bc`; echo $a
{	mkdir X$i && cd X$i
	pp=`echo "scale=4;$i/125.47"|bc`
	java -Xmx12g -jar $picard DownsampleSam I=$sIn O=X$i.SZ-1.filtered.bam P=$pp 
##ln -s ../mhk.dedup.bam.bai ./0.05_mhk.dedup.bam.bai #The old index file can't be used and need to create a new.
	samtools index X$i.SZ-1.filtered.bam
	$rohan --rohmu 2e-5 -t 8 --size 50000 --step 100 -o X$i $ref $sBam_path/X$i/X$i.SZ-1.filtered.bam
	cd ..
}&
done
wait



