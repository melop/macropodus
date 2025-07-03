for sD in ./*.hEst.gz; do
{	sF=$(basename $sD)
	sP=${sF/.hEst.gz/}
	echo $sP
	less $sD | awk -v OFS="\t" '{if(substr($1,1,1)!="#"){print $1,$2,$3,$5}}' > $sP.rohan.het.bed
}&
done
wait
