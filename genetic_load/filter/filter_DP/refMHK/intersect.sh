bedtools intersect -a allpops.gvcf.depth.dedup.bed -b refMHK.cds.bed -wa > allpops.cds.depth.bed
sed -i 's/\./0/g' allpops.cds.depth.bed

#有些基因存在重叠区间导致结果中某些位点出现重复记录，因此后续还需要dedup.R去除重复
