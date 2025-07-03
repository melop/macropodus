vcfcaller=/data/software/msmc-tools/vcfAllSiteParser.py

CONVERT=/data/software/msmc-tools/generate_multihetsep.py
BCFTOOLS="bcftools"
sVCFFolder=/data2/projects/zwang/m.hk/all_gvcf
sVCFSuffix=.genotyped.g.vcf.gz
sChrSuffix="scf_"

sRef=mhk
arrSamples=( MOPld-5 MOPdxs-1 GL-4 MOPhn-1 YN-mop1 HK-1 MHKmls-1 MHKmls-6 JX-4.bam )
arrCovUpper=( 66 104 13 117 23 18 58 73 31 )
arrCovLower=( 5 5 5 5 5 5 5 5 5 )

arrChr=( `seq 1 23` )

sOutDir=formsmc2_in/


for nSample in "${!arrSamples[@]}"; do 
	sSample=${arrSamples[$nSample]};
	nCovUpper=${arrCovUpper[$nSample]};
	nCovLower=${arrCovLower[$nSample]};
	sFilterCmd=${arrFilterCmd[$nSample]};

	for nChr in "${arrChr[@]}"; do 
		sChr=$sRef$sChrSuffix$nChr
		sChrDir=$sOutDir/$sSample/$sChr
		mkdir -p $sChrDir
		( if [ ! -s $sChrDir/out_mask.bed.gz ];then $BCFTOOLS view -r "$sChr" -e 'INFO/DP>'$nCovUpper' || INFO/DP<'$nCovLower'' $sVCFFolder/$sSample$sVCFSuffix | $vcfcaller $sChr $sChrDir/out_mask.bed.gz | gzip -c > $sChrDir/out.vcf.gz; fi; \
		if [ ! -e $sChrDir/done.txt ]; then $CONVERT --mask $sChrDir/out_mask.bed.gz $sChrDir/out.vcf.gz > $sChrDir/formsmc2.multihetsep.txt 2> $sChrDir/formsmc2.multihetsep.log && touch $sChrDir/done.txt; fi; ) 2>$sChrDir/log.txt &
	done

done

wait

