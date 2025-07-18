#SV_MHK=/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/SV.w.indelmodifier.polarized.out.txt
#SV_MOP=/data2/projects/zwang/m.op/Snpeff/all_MHK_MOP/rho_V2.0/SV.w.indelmodifier.polarized.out.txt
SV_MHK=refMHK.dedup.SV.w.indelmodifier.polarized.out.txt
SV_MOP=refMOP.dedup.SV.w.indelmodifier.polarized.out.txt
sPop=allpops

#step1: stats AF with statsAF.R

#step2:
#awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$10,$11,$14}' /data2/projects/zwang/macropodus_compare/genetic_load/all_MHK_MOP/liftover/$sPop.cds.bothDP.bed > temp.$sPop.cds.bothDP.bed
#bedtools intersect -a temp.$sPop.cds.bothDP.bed -b $sPop.refMOP.SV.withAF.bed -wa -wb > $sPop.cds.bothDP.MOPAF.bed
#awk -v OFS='\t' '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14}' $sPop.cds.bothDP.MOPAF.bed > temp.$sPop.cds.bothDP.MOPAF.bed
#bedtools intersect -a temp.$sPop.cds.bothDP.MOPAF.bed -b $sPop.refMHK.SV.withAF.bed -wa -wb> $sPop.cds.bothDP.bothAF.bed

#step3:filter DP and SNP with filter.R

#step4:
awk -F"\t" -v OFS='\t' '{ start=$2-1;  printf("%s",$1"\t"start"\t"$2"\t");$1="";$2="";print}' $SV_MHK > $sPop.refMHK.SV.polarized.out.bed
awk -F"\t" -v OFS='\t' '{ start=$2-1;  printf("%s",$1"\t"start"\t"$2"\t");$1="";$2="";print}' $SV_MOP > $sPop.refMOP.SV.polarized.out.bed
bedtools intersect -a $sPop.refMHK.SV.polarized.out.bed -b $sPop.refMHK.SV.filtered.bed > $sPop.refMHK.SV.polarized.out.filtered.bed
bedtools intersect -a $sPop.refMOP.SV.polarized.out.bed -b $sPop.refMOP.SV.filtered.bed > $sPop.refMOP.SV.polarized.out.filtered.bed
