FILE=../../MHK.noMissing.treemix.frq.gz
treemix=/data/software/treemix/src/treemix
Migrate=1

$treemix -i $FILE -k 1000 -root MHKhn -m $Migrate -tf ../ML/reroot.txt -se  -o migrate_${Migrate} > log.txt 2>&1

