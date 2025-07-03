for i in ../ensembl/*; do
	if [ -d $i ]; then
		sp=`basename $i`
		zcat $i/*.pep.all.fa.gz | php /data/software/killigenomics/annotation/UPhO_Ensembl/preprocessAA.php > $i.filtered.fa 2> $i.filtered.tab
		/data/software/UPhO/minreID.py $i.filtered.fa $sp \|
	fi
done

cat *.fst > allENSproteins.fa
