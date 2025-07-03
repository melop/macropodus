sD=/data2/projects/zwang/m.hk/ROH/ROHan/runs_of_het
sD2=/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2/seperate
gff=/data/projects/rcui/mhk/annotations/funannotate/train/mhk.longest_isoform.genesymbol.gff3

##balancing signature in shared hetrozygosity of three high Froh populations
#bedtools intersect -a $sD/MHKjx/MHKjx.het.union.bed -b $sD/MHKmlh_qns/MHKmlh_qns.het.union.bed > temp.bed
#bedtools intersect -a temp.bed -b $sD/MOPfq_pt/MOPfq_pt.het.union.bed > MHKjx.MHKmlh_qns.MOPfq_pt.intersected.HET.bed

#bedtools intersect -a MHKjx.MHKmlh_qns.MOPfq_pt.intersected.HET.bed -b $sD2/MHKjx.balanced.bed > temp2.bed
#bedtools intersect -a temp2.bed -b $sD2/MHKmlh_qns.balanced.bed > temp3.bed
#bedtools intersect -a temp3.bed -b $sD2/MOPfq_pt.balanced.bed > MHKjx.MHKmlh_qns.MOPfq_pt.intersected.HET.balanced.bed

#balancing signature in shared hetrozygosity of two high Froh populations
list=("MHKjx" "MHKmlh_qns" "MOPfq_pt")
for (( i=0; i<${#list[@]}; i++ ));do
	for (( j=$i+1; j<${#list[@]}; j++ ));do
		(bedtools intersect -a $sD/${list[i]}/${list[i]}.het.union.bed -b $sD/${list[j]}/${list[j]}.het.union.bed > ${list[i]}.${list[j]}.intersected.HET.bed
		bedtools intersect -a ${list[i]}.${list[j]}.intersected.HET.bed -b $sD2/${list[i]}.balanced.bed >  ${list[i]}.${list[j]}.temp.bed
		bedtools intersect -a ${list[i]}.${list[j]}.temp.bed -b $sD2/${list[j]}.balanced.bed > ${list[i]}.${list[j]}.intersected.HET.balanced.bed
		while read line; do
			chrom=$(echo $line | awk '{print $1}')
			start=$(echo $line | awk '{print $2+1}')
			end=$(echo $line | awk '{print $3}')
			awk -v OFS="\t" -v chrom=$chrom -v start=$start -v end=$end '{if($1==chrom && $3=="gene" && (($4 > start-1 && $4 < end+1) || ($5 > start-1 && $5 < end+1))){print $0}}' $gff >> ${list[i]}.${list[j]}.intersected.HET.balanced.gff3
		done < ${list[i]}.${list[j]}.intersected.HET.balanced.bed )&
	done
done
wait
