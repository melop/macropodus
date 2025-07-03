awk -v OFS="\t" '{print $1,$2}' pairedMHKMOP.dxy.top0.01.txt | sort | uniq > paired.list.txt
cat paired.list.txt | while read line
do
    a=$(echo $line | awk '{print $1}');
    b=$(echo $line | awk '{print $2}');
    touch $a.$b.dxy.bed
    awk -v OFS="\t" '$1 == "'$a'" && $2 == "'$b'" {print $3,$4-1,$5,$6}' pairedMHKMOP.dxy.top0.01.txt > $a.$b.dxy.bed
done

