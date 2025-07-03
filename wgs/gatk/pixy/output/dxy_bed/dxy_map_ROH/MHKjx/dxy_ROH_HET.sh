module load bedtools
comb=MOPdsy.MHKjx
for i in JX-1 JX-2 JX-3 JX-4;do
	het=$i.het.bed
	hom=$i.hom.bed
#	echo $het
#	echo $hom
	bedtools intersect -a $comb.dxy.bed -b $het -wa > $i.het.dxy.bed
	bedtools intersect -a $comb.dxy.bed -b $hom -wa > $i.hom.dxy.bed
	awk -v OFS="\t" 'BEGIN{row=0;dxy=0};{row+=1;dxy+=$4};END{print "'$i'","meandxyonHET",dxy/row}' $i.het.dxy.bed >> $comb.meandxy.HET.HOM.txt
	awk -v OFS="\t" 'BEGIN{row=0;dxy=0};{row+=1;dxy+=$4};END{print "'$i'","meandxyonHOM",dxy/row}' $i.hom.dxy.bed >> $comb.meandxy.HET.HOM.txt
done
