#mkfifo _tab.csv
#cat BLAST_results_allproteins.csv | tr "," "\t" > _tab.csv &
python2 /data/software/UPhO/BlastResultsCluster.py -in BLAST_results_allproteins.csv  -mcl -e 1e-10
