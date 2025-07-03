#cp /data/projects/zwang/ROHan_test/0.05_SZ1_on_mhk.h.txt ./ # this is rohan estimated heterozygosity at low coverage
#cat 0.05_SZ1_on_mhk.h.txt | cut -d" " -f1,2,3,4 | grep -v '#' | awk '{print $1"\t"$2"\t"$3"\t"$4}' > rohan.het.bed
#cat ../mhk.hom.win.tsv | awk '{if ($1!="chr") {print $1"\t"$3"\t"$4} }' > roh.bed
#bedtools intersect -a rohan.het.bed -b roh.bed -wa -f 0.99 > roh.wins.bed


cat ./mhk.hom.win.tsv | awk '{if ($1!="chr") {print $1"\t"$3"\t"$4} }' > roh.bed
for sP in /data2/projects/zwang/m.hk/ROH/simulation1.2/X*; do
{	sD=$(basename $sP)
#	echo $sD
	less $sP/*.hEst.gz | awk -v OFS="\t" '{if(substr($1,1,1)!="#"){print $1,$2,$3,$5}}' > $sD.rohan.het.bed
	bedtools intersect -a $sD.rohan.het.bed -b roh.bed -wa -f 0.99 > $sD.roh.wins.bed
}&
done
wait
