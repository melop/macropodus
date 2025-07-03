FILE=../../MHK.noMissing.treemix.frq.gz
treemix=/data/software/treemix/src/treemix

$treemix -i $FILE -k 1000 --root MHKhn -o mltree > log.txt 2>&1
