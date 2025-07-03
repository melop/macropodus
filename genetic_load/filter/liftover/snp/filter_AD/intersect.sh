bedtools intersect -a allpops.refMHK.gvcf.AD.bed -b temp.refMHK.bed -wa > allpops.refMHK.cds.AD.bed &
bedtools intersect -a allpops.refMOP.gvcf.AD.bed -b temp.refMOP.bed -wa > allpops.refMOP.cds.AD.bed &
wait
