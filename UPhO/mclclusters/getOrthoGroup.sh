MINTAXA=3
MINBOOSTRAP=0.50
PROT=allproteins.fa

mkdir -p trees
for i in UPhO_Seqs/*bipartition*.out; do
  sGroup=`echo "$i" | grep -oP "Group_\d+"`
  echo $sGroup
  cp $i trees/$sGroup.tre
done

rm UPhO_orthogroups.csv
rm -R UPhO_branches/

python2 /data/software/UPhO/UPhO.py  -in trees/*.tre -m $MINTAXA -S $MINBOOSTRAP -ouT -iP -R ../$PROT > ortholog.log 2>&1;

