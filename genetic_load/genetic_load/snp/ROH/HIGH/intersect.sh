sBed=/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/ROH

for i in $sBed/MHKqns* $sBed/SZ* $sBed/FJ* $sBed/PT* $sBed/JX*; do
{	sF=$(basename $i)
	ID=${sF/.bed/}
#	echo $ID
	bedtools intersect -a allpops.refMHK.NONSYN.polarized.out.filtered.withAD.dedup.withConsurf.dedup.bed -b $i > allpops.refMHK.NONSYN.filtered.$ID.bed 
} &
done
wait
