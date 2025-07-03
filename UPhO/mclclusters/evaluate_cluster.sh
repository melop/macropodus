echo "" > evaluate_cluster.out.txt

for i in out.seq.mci.I*; do
	php evaluate_cluster.php $i 18 >> evaluate_cluster.out.txt
done
