liftOver=/data/software/LiftOver/liftOver

cp /data2/projects/zwang/macropodus_compare/synteny/mop2mhk/mop2mhk.chain ./renamed.mop2mhk.chain 
sed -i "s/mhkscf_/chr/g" renamed.mop2mhk.chain
sed -i "s/mopscf_/chr/g" renamed.mop2mhk.chain

#pop=allpops
sF=/data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/filter_DP/refMOP/allpops.cds.depth.dedup.bed
awk -v OFS="\t" '{print $0,$1,$2,$3}' $sF > temp.allpops_on_MOP.cds.depth.bed
sed -i "s/^mopscf_/chr/" temp.allpops_on_MOP.cds.depth.bed

$liftOver temp.allpops_on_MOP.cds.depth.bed renamed.mop2mhk.chain allpops_on_MOP2MHK.withOriginName.cds.depth.bed unmapped.txt
