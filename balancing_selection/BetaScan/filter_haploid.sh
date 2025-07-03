sPop=MOPdsy
sPath=/data2/projects/zwang/m.hk/ROH/BetaScan/$sPop/beta_input_for_chrom
for sP in $sPath/*.beta.txt; do
{	sF=$(basename $sP)
	sStem=${sF/.beta.txt/}
	echo $sStem
	less $sP | awk -v OFS="\t" -v sStem=$sStem 'BEGIN{left=0;right=0};{if($1-250<0){left=0}else{left=$1-250};right=$1+250;print(sStem,left,right,$2,$3)}' > $sStem.beta.bed
	less $sP | awk -v OFS="\t" -v sStem=$sStem 'BEGIN{left=0;right=0};{if($3<3){if($1-250<0){left=0}else{left=$1-250};right=$1+250;print(sStem,left,right)}}' > $sStem.dis.beta.bed
	bedtools subtract -A -a $sStem.beta.bed -b $sStem.dis.beta.bed > $sStem.filtered.beta.bed
	less $sStem.filtered.beta.bed | awk -v OFS="\t" '{Origin=$3-250;print(Origin,$4,$5)}' > $sStem.filtered.beta.txt
#	rm $sStem.beta.bed $sStem.dis.beta.bed $sStem.filtered.beta.bed
} &
done
wait

