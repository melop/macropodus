for i in ./*.betascores.txt; do
	less $i | awk '{if(NR>1){print $2}}' >> temp.all_chrom.betascores.txt
done

awk '{ print $0 }' temp.all_chrom.betascores.txt | sort -n > all_chrom.betascores.txt


awk '{ a[NR] = $1 };END{ n = int(NR*0.999)+1;print "top 0.1 percentile",a[n]}' all_chrom.betascores.txt
awk '{ a[NR] = $1 };END{ n = int(NR*0.995)+1;print "top 0.5 percentile",a[n]}' all_chrom.betascores.txt
awk '{ a[NR] = $1 };END{ n = int(NR*0.99)+1;print "top 1 percentile",a[n]}' all_chrom.betascores.txt
