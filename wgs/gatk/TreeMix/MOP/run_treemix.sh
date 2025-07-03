treemix=/data/software/treemix/src/treemix
plink=/data/software/Plink/plink
Na=MOP #Ref 
##step1:remove sites with missing data because treemix does not like missing data
#vcftools --gzvcf ./all.mop_selected_refmhk_filtersnps.vcf.gz --max-missing 1 --recode --stdout | gzip > $Na.noMissing.vcf.gz;

Na=MOP.noMissing
#step2:create clust file
#bcftools query -l ./all.mop_selected_refmhk_filtersnps.vcf.gz | awk '{split($1,pop,"."); print $1"\t"$1"\t"pop[2]}' > $Na.cluster; #Check and change this file
#step3:transfer to format for treemix using scripts in treemix
bash /data2/projects/zwang/m.hk/Geneflow/TreeMix/MOP/vcf2treemix.sh $Na.vcf.gz $Na.cluster
#step4:assumes 0-8 geneflow and run treemix
mkdir out
cd out
for m in {0..10}; do
	for i in {1..6}; do    
	(	$treemix -se -bootstrap \
	    	-i ../$Na.treemix.frq.gz \
	    	-o $Na.${m}.${i} \
	    	-root MOPyn \
	    	-m ${m} \
	    	-k 1000 \
	    	-noss > treemix.$m.$i.log
	) &
	done
done
wait

