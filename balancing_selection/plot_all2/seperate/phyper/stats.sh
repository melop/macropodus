sD=/data2/projects/zwang/m.hk/ROH/BetaScan/plot_all2
sD2=/data2/projects/zwang/m.hk/ROH/ROHan/runs_of_het
list=("MHKhn" "MHKjx" "MHKmlh_qns" "MHKmls_pop1" "MHKmls_pop2" "MOPdsy" "MOPdxs_sg" "MOPfq_pt")
#for f in $sD/*.beta.bed; do
#	sF=$(basename $f)
#	sStem=${sF/.beta.bed/}
#	awk -v OFS="\t" '{if($4>0){print $0}}' $f > $sStem.balanced.bed
#done


for (( i=0; i<${#list[@]}; i++ ));do
	for (( j=$i; j<${#list[@]}; j++ ));do
#		echo "${list[i]} ${list[j]}"
		{
#		awk -v OFS="\t" '{if($4>0){print $0}}' $sD/${list[i]}.beta.bed > ${list[i]}.balanced.bed
#		awk -v OFS="\t" '{if($4>0){print $0}}' $sD/${list[j]}.beta.bed > ${list[j]}.balanced.bed
#		bedtools intersect -a $sD2/${list[i]}/${list[i]}.het.union.bed -b $sD2/${list[j]}/${list[j]}.het.union.bed > ${list[i]}.${list[j]}.intersected.HET.bed

		bedtools intersect -a ../${list[i]}.${list[j]}.intersected.HET.bed -b ../${list[i]}.balanced.bed -wa | sort -k1,1 -k2,2n -k3,3n | uniq > ${list[i]}.${list[j]}.${list[i]}.balanced.bed #SP1存在平衡选择的共有杂合区间
		bedtools intersect -a ../${list[i]}.${list[j]}.intersected.HET.bed -b ../${list[j]}.balanced.bed -wa | sort -k1,1 -k2,2n -k3,3n | uniq > ${list[i]}.${list[j]}.${list[j]}.balanced.bed #SP2存在平衡选择的共有杂合区间

		bedtools intersect -a ${list[i]}.${list[j]}.${list[i]}.balanced.bed -b ../${list[j]}.balanced.bed -wa | sort -k1,1 -k2,2n -k3,3n | uniq > ${list[i]}.${list[j]}.${list[i]}.${list[j]}.bothbalanced.bed #SP1，SP2同时存在平衡选择的共有杂合区间
		bedtools subtract -a ${list[i]}.${list[j]}.${list[i]}.balanced.bed -b ../${list[j]}.balanced.bed -A > ${list[i]}.${list[j]}.only.${list[i]}.balanced.bed #仅SP1存在平衡选择的共有杂合区间
		bedtools subtract -a ${list[i]}.${list[j]}.${list[j]}.balanced.bed -b ../${list[i]}.balanced.bed -A > ${list[i]}.${list[j]}.only.${list[j]}.balanced.bed #仅SP2存在平衡选择的共有杂合区间

		bedtools subtract -a ../${list[i]}.${list[j]}.intersected.HET.bed -b ../${list[i]}.balanced.bed -A > ${list[i]}.${list[j]}.no.${list[i]}.balanced.bed #SP1不存在平衡选择的共有杂合区间
		bedtools subtract -a ${list[i]}.${list[j]}.no.${list[i]}.balanced.bed -b ../${list[j]}.balanced.bed -A > ${list[i]}.${list[j]}.${list[i]}.${list[j]}.bothno.balanced.bed #SP1,SP2均不存在平衡选择的共有杂合区间
		}&
	done
done
wait
