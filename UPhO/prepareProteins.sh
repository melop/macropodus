for i in ./macropodus_proteins/*; do
	if [ -d $i ]; then
		sp=`basename $i`
		zcat -f $i/*.fa | php preprocessAA.php > $i.filtered.fa 2> $i.filtered.tab
		/data/software/UPhO/minreID.py $i.filtered.fa $sp \|
	fi
done

cat *.fst > allmacropodusproteins.fa
