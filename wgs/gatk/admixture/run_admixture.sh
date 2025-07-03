Name=all_refmhk_new
admixture=/data/software/admixture_linux-1.3.0/admixture
for K in 3 4 5 6 7 8 9 10 11 12 13 14 15; do
	$admixture --cv $Name.bed $K -j48 | tee log${K}.out
done
