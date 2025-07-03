bed=/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/liftover/snp/allpops.cds.bothDP.without_discarded.bed
MOP=allpops.refMOP.cds.AD.filtered.bed
MHK=allpops.refMHK.cds.AD.filtered.bed

bedtools intersect -a $bed -b $MOP -wa > temp1.allpops.cds.bothDP.without_discarded.filterAD.bed

awk -v OFS="\t" '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14}' temp1.allpops.cds.bothDP.without_discarded.filterAD.bed > temp2.allpops.cds.bothDP.without_discarded.filterAD.bed

bedtools intersect -a temp2.allpops.cds.bothDP.without_discarded.filterAD.bed -b $MHK -wa > temp3.allpops.cds.bothDP.without_discarded.filterAD.bed

awk -v OFS="\t" '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14}' temp3.allpops.cds.bothDP.without_discarded.filterAD.bed > allpops.cds.bothDP.without_discarded.filterAD.bed

rm temp1* temp2* temp3*
