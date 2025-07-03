awk -v OFS="\t" '!seen[$1, $2, $3]++' allpops.gvcf.depth.bed > allpops.gvcf.depth.dedup.bed
