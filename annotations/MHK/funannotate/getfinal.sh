sp=mhk
genome=trained/predict_misc/genome.softmasked.fa

grep -vP "^#" trained/update_misc/pasa/Macropodus_hongkongensis_SZDP2020_pasa.gene_structures_post_PASA_updates.1835267.addedback.gff3 | \
sed 's/FUN_/'$sp'/g' | \
sed 's/novel_/'$sp'_/g' \
> $sp.full.gff3



grep -P "\ttRNA\t" trained/update_results/Macropodus_hongkongensis_SZDP2020.gff3 > $sp.tRNA.gff3

/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_single_longest_isoform.pl $sp.full.gff3 > $sp.longest_isoform.gff3

/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl $sp.full.gff3 $genome cDNA > $sp.full.mrna.fa
/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl $sp.full.gff3 $genome prot > $sp.full.prot.fa
/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl $sp.full.gff3 $genome CDS > $sp.full.cds.fa

/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl $sp.longest_isoform.gff3 $genome cDNA > $sp.longest_isoform.mrna.fa
/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl $sp.longest_isoform.gff3 $genome prot > $sp.longest_isoform.prot.fa
/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl $sp.longest_isoform.gff3 $genome CDS > $sp.longest_isoform.cds.fa

