less DSY-1.hEst.gz | awk -v OFS="\t" '{if(substr($1,1,1)!="#"){print $1,$2,$3,$5}}' > DSY-1.rohan.het.bed
