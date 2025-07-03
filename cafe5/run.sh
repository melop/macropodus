
CAFE=/data/software/CAFE5/bin/cafe5
CPU=96

for i in `seq 6 10`; do
	mkdir -p runresults/$i
$CAFE -c $CPU -t tree.tre -i orthogroups.genecount.txt -p -eresults/Base_error_model.txt -k 3 -o  runresults/$i >  runresults/$i/est.log 2>&1 
done
