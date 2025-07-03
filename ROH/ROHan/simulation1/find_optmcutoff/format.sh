sPath="/data/projects/zwang/ROHan_test"
for sF in $sPath/*.gz; do
	sH=$(basename $sF)
	sStem=${sH/.hEst.gz/}
#	echo $sStem
	less $sF | awk '{print $1,$2,$3,$5}' > $sStem.h.txt
done
#less SZ-1.h.ZWANG.txt | awk 'NR==1 {print $0}' > SZ-1.h.RC.txt
#less $sIn2 | awk '{print $1,$2,$3,$9}' >> SZ-1.h.RC.txt

