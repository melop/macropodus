BetaScan=/data/software/BetaScan/BetaScan.py
sPop=MOPdsy
sPath=/data2/projects/zwang/m.hk/ROH/BetaScan/$sPop/beta_input_for_chrom

for i in $sPath/*.filtered.beta.txt; do
{	sF=$(basename $i)
	sStem=${sF/.filtered.beta.txt/}
	echo $sStem
	python $BetaScan -i $i -fold -o $sPop.$sStem.betascores.txt
} &
done
wait
