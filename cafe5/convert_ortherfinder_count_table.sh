IN=/data2/projects/zwang/macropodus_compare/synteny/genespace_RC/rundir/orthofinder/Results_Nov01/Orthogroups/Orthogroups.GeneCount.tsv
cat $IN | cut -f1-11 | awk '{if ($1=="Orthogroup") {print "Desc\t"$0 } else {print "(null)\t"$0 } }' > orthogroups.genecount.txt
