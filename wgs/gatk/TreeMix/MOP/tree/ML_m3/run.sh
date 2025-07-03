FILE=../../MOP.noMissing.treemix.frq.gz
treemix=/data/software/treemix/src/treemix
Migrate=3

$treemix -i $FILE -k 1000 -root MOPyn -m $Migrate -tf ../ML/reroot.txt -se  -o migrate_${Migrate} > log.txt 2>&1

