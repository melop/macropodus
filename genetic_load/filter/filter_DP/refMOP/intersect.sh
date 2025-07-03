bedtools intersect -a allpops.gvcf.depth.dedup.bed -b refMOP.cds.bed -wa > allpops.cds.depth.bed
sed -i 's/\./0/g' allpops.cds.depth.bed
