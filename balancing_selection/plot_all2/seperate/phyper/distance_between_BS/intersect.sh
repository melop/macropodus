list=("MHKhn" "MHKjx" "MHKmlh_qns" "MHKmls_pop1" "MHKmls_pop2" "MOPdsy" "MOPdxs_sg" "MOPfq_pt")

for (( i=0; i<${#list[@]}; i++ ));do
        for (( j=$i; j<${#list[@]}; j++ ));do
	{	bedtools intersect -a ${list[i]}.BS.bed -b ../${list[i]}.${list[j]}.${list[i]}.${list[j]}.bothbalanced.bed -wa -wb > ${list[i]}.BS.${list[i]}.${list[j]}.shearedHET.bed
		bedtools intersect -a ${list[j]}.BS.bed -b ../${list[i]}.${list[j]}.${list[i]}.${list[j]}.bothbalanced.bed -wa -wb > ${list[j]}.BS.${list[i]}.${list[j]}.shearedHET.bed
	} &
	done
done
wait
