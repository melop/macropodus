sD="/data2/projects/zwang/m.hk/all_gvcf"
echo -e "MHK\t/data/projects/zwang/m.hk/GATK/mhk.LG.nextdenovo.sgs.2.sorted.fa" > ./alltaxa.txt
for i in $sD/MOPhn-6*gz $sD/MHKmls-6*gz $sD/MOPld-4*gz $sD/MHKqns-2*gz $sD/MOPsg-2*gz $sD/DSY-1*gz $sD/YN-mop1*gz; do
	sSample=$(basename $i)
	sStem=${sSample%%.*}
	echo -e "$sStem\t$i" >> ./alltaxa.txt
done
