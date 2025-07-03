mkdir macropodus_proteins
cd macropodus_proteins

mkdir MHK
cat /data/projects/rcui//mhk/annotations/funannotate/train/mhk.longest_isoform.prot.fa | sed 's/*//g' > ./MHK/prot.fa

mkdir MOP
cat /data/projects/rcui//mop/annotations/funannotate/train/mop.longest_isoform.prot.fa | sed 's/*//g' > ./MOP/prot.fa
