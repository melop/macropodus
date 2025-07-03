sP=/data2/projects/zwang/m.hk/ROH/ROHan/runs_of_hom

for sD in $sP/*/*ROH.bed; do
	sF=$(basename $sD)
	sStem=${sF%%.*}
#	echo $sStem
	LR=429175390 #用于统计的参考基因组长度
	awk -v LR="$LR" 'BEGIN{ROH=0};{if($3-$2>50000){ROH+=$3-$2}};END{printf("'$sStem' Froh: %10.5f\n",ROH/LR)}' $sD >> allpops.Froh_100kb.txt
done

