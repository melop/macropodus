for sP in /data2/projects/zwang/m.hk/ROH/sw_rohan_DSY1/simulation1.2/X*; do
{       sD=$(basename $sP)
#       echo $sD
        less $sP/*.hEst.gz | awk -v OFS="\t" '\
                {if(substr($1,1,6) == "mhkscf")\
                        {split($1,arr,"_");if (arr[2]<24 && $4>=25000 && $5<=0.00082){print $1,$2,$3,$5}}}' >> $sD.hom.bed
        bedtools merge -i $sD.hom.bed > $sD.ROH.bed
	rm $sD.hom.bed
}&
done
wait

