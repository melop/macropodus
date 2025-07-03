sD=/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2
sD2=/data2/projects/zwang/m.hk/ROH/ROHan/runs_of_het
list=("MHKhn" "MHKjx" "MHKmlh_qns" "MHKmls_pop1" "MHKmls_pop2" "MOPdsy" "MOPdxs_sg" "MOPfq_pt")
for f in $sD/*.beta.bed; do
	sF=$(basename $f)
	sStem=${sF/.beta.bed/}
	awk -v OFS="\t" '{if($4>0){print $0}}' $f > $sStem.balanced.bed
done


for (( i=0; i<${#list[@]}; i++ ));do
	for (( j=$i; j<${#list[@]}; j++ ));do
#		echo "${list[i]} ${list[j]}"
		{
#		awk -v OFS="\t" '{if($4>0){print $0}}' $sD/${list[i]}.beta.bed > ${list[i]}.balanced.bed
#		awk -v OFS="\t" '{if($4>0){print $0}}' $sD/${list[j]}.beta.bed > ${list[j]}.balanced.bed
		bedtools intersect -a $sD2/${list[i]}/${list[i]}.het.union.bed -b $sD2/${list[j]}/${list[j]}.het.union.bed > ${list[i]}.${list[j]}.intersected.HET.bed

		bedtools intersect -a ${list[i]}.${list[j]}.intersected.HET.bed -b ${list[i]}.balanced.bed -wa | sort -k1,1 -k2,2n -k3,3n | uniq > ${list[i]}.${list[j]}.${list[i]}.balanced.bed
		bedtools intersect -a ${list[i]}.${list[j]}.${list[i]}.balanced.bed -b ${list[j]}.balanced.bed -wa | sort -k1,1 -k2,2n -k3,3n | uniq > ${list[i]}.${list[j]}.${list[i]}.${list[j]}.bothbalanced.bed

		bedtools subtract -a ${list[i]}.${list[j]}.${list[i]}.${list[j]}.bothbalanced.bed -b ${list[i]}.balanced.bed | awk -v OFS="\t" '{print $0,0}' >  ${list[i]}.${list[j]}.${list[i]}.subtracted.bed
		bedtools subtract -a ${list[i]}.${list[j]}.${list[i]}.${list[j]}.bothbalanced.bed -b ${list[j]}.balanced.bed | awk -v OFS="\t" '{print $0,0}' >  ${list[i]}.${list[j]}.${list[j]}.subtracted.bed

		cat ${list[i]}.balanced.bed ${list[i]}.${list[j]}.${list[i]}.subtracted.bed | sort -k1,1 -k2,2n > ${list[i]}.${list[j]}.${list[i]}.seperated.bed
		cat ${list[j]}.balanced.bed ${list[i]}.${list[j]}.${list[j]}.subtracted.bed | sort -k1,1 -k2,2n > ${list[i]}.${list[j]}.${list[j]}.seperated.bed
		}&
	done
done
wait
