module load bedtools
comb=MOPdsy.MHKmlh_qns
for i in SZ-1 SZ-2 SZ-3 SZ-4 MHKqns-1 MHKqns-2 MHKqns-3;do
	het=$i.het.bed
	hom=$i.hom.bed
#	echo $het
#	echo $hom
	bedtools intersect -a $comb.dxy.bed -b $het -wa > $i.het.dxy.bed
	bedtools intersect -a $comb.dxy.bed -b $hom -wa > $i.hom.dxy.bed
done
