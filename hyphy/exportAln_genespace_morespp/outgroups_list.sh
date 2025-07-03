sPath1="/data/projects/zwang/ncbi"
sPath2="/data2/projects/zwang/ncbi"
echo -e "#name\tcds" > outgroups.txt
for i in $sPath1/* $sPath2/*; do
	if [ -d $i ]; then
		sStem=$(basename $i)
		echo -e "$sStem\t$i/longest_isoform.cds.fa" >> outgroups.txt
	fi
done

