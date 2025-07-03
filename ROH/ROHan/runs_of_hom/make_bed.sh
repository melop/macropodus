sPath=/data2/projects/zwang/m.hk/ROH/ROHan/all_het

for i in $sPath/FJ-1 $sPath/PT-1 $sPath/PT-2; do
	sH=$i/*.hEst.gz
#	echo $sH
	sStem=$(basename $i)
#	echo $sStem
#	echo -e "#CHROM\tBEGIN\tEND\th" > $sStem.het.bed
	less $sH | awk -v OFS="\t" '\
		{if(substr($1,1,6) == "mhkscf")\
			{split($1,arr,"_");if (arr[2]<24 && $4>=25000 && $5<=0.00082){print $1,$2,$3,$5}}}' >> $sStem.hom.bed
	bedtools merge -i $sStem.hom.bed > $sStem.ROH.bed
done

cat *.ROH.bed > MOPfq_pt.ROH.bed
