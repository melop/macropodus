sBed=/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/ROH

for i in $sBed/FJ*.HET.bed $sBed/PT*.HET.bed; do
	sF1=$(basename $i)
	ID1=${sF1/.HET.bed/}
	for j in $sBed/MHKqns*.HET.bed $sBed/SZ*.HET.bed $sBed/JX*.HET.bed; do
		sF2=$(basename $j)
		ID2=${sF2/.HET.bed/}
		bedtools intersect -a $i -b $j > $ID1.$ID2.HET.bed
		bedtools intersect -a allpops.refMHK.NONSYN.polarized.out.filtered.withAD.dedup.withConsurf.dedup.bed -b $ID1.$ID2.HET.bed > allpops.refMHK.NONSYN.filtered.$ID1.$ID2.HET.bed
	done
done

for i in $sBed/FJ*.ROH.bed $sBed/PT*.ROH.bed; do
        sF1=$(basename $i)
        ID1=${sF1/.ROH.bed/}
        for j in $sBed/MHKqns*.ROH.bed $sBed/SZ*.ROH.bed $sBed/JX*.ROH.bed; do
                sF2=$(basename $j)
                ID2=${sF2/.ROH.bed/}
                bedtools intersect -a $i -b $j > $ID1.$ID2.ROH.bed
                bedtools intersect -a allpops.refMHK.NONSYN.polarized.out.filtered.withAD.dedup.withConsurf.dedup.bed -b $ID1.$ID2.ROH.bed > allpops.refMHK.NONSYN.filtered.$ID1.$ID2.ROH.bed
        done
done

#for i in $sBed/MHKqns* $sBed/SZ* $sBed/FJ* $sBed/PT* $sBed/JX*; do
#{	sF=$(basename $i)
#	ID=${sF/.bed/}
##	echo $ID
#	bedtools intersect -a allpops.refMHK.SV.polarized.out.filtered.bed -b $i > allpops.refMHK.SV.filtered.$ID.bed 
#} &
#done
#wait
