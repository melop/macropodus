sD=/data2/projects/zwang/m.hk/ROH/ROHan/runs_of_het
list=("MHKhn" "MHKjx" "MHKmlh_qns" "MHKmls_pop1" "MHKmls_pop2" "MOPdsy" "MOPdxs_sg" "MOPfq_pt")
for (( i=0; i<${#list[@]}; i++ ));do
	for (( j=$i; j<${#list[@]}; j++ ));do
#		echo "${list[i]} ${list[j]}"
	(	bedtools intersect -a $sD/${list[i]}/${list[i]}.het.union.bed -b $sD/${list[j]}/${list[j]}.het.union.bed > ${list[i]}.${list[j]}.intersected.HET.bed

		# 找到文件A与文件B的交集，并求和第四列
		bedtools intersect -a ${list[i]}.${list[j]}.intersected.HET.bed -b ${list[i]}.beta.bed -wa -wb | \
			          awk -v OFS='\t' '{sum[$1":"$2"-"$3]+=$7} END {for (i in sum) print i, sum[i]}' | \
				              sort -k1,1V -k2,2n > ${list[i]}.${list[j]}.${list[i]}.singlebeta.bed

		# 找到文件A与文件C的交集，并求和第四列
		bedtools intersect -a ${list[i]}.${list[j]}.intersected.HET.bed -b ${list[j]}.beta.bed -wa -wb | \
			          awk -v OFS='\t' '{sum[$1":"$2"-"$3]+=$7} END {for (i in sum) print i, sum[i]}' | \
				              sort -k1,1V -k2,2n > ${list[i]}.${list[j]}.${list[j]}.singlebeta.bed

		# 合并文件A与文件B、文件C的交集，并添加第四列和第五列求和结果
		join -a 1 -a 2 -e 0 -o 0,1.2,2.2 ${list[i]}.${list[j]}.${list[i]}.singlebeta.bed ${list[i]}.${list[j]}.${list[j]}.singlebeta.bed > ${list[i]}.${list[j]}.${list[i]}.${list[j]}.bothbeta.bed	)&
	done
done
wait
