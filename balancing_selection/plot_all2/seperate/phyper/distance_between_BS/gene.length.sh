less mhk.longest_isoform.genesymbol.gff3 | awk -v OFS="\t" '{if($3=="gene"){print $1,$4,$5,$5-$4}}' > gff.genes.length.txt
