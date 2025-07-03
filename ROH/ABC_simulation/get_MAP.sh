# cat out/* | awk '{if ($4>=5.25 && $4<=5.32 && $5>=0.25 && $5<=0.35 && $6>0.72 && $6<0.83) print}' | awk '{gen+=$1; pop+=$2; rec+=$3} END{print NR, gen/NR, pop/NR, rec/NR}'

cat out/* | head -n1 > abc.reps.txt
cat out/* | grep -v gen >>  abc.reps.txt
