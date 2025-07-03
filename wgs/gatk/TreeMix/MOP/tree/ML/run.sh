FILE=../../MOP.noMissing.treemix.frq.gz
treemix=/data/software/treemix/src/treemix

$treemix -i $FILE -k 1000 --root MOPyn -o mltree > log.txt 2>&1
