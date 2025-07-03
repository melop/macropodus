outstem=joincalled
genome=/public3/group_crf/home/g20wangzhx36/m.hk/nextdenovo/mhk.LG.nextdenovo.sgs.2.sorted.fa
gatk=/public3/group_crf/home/wanglinlong8405/genome/lty_genome/resequencing/01.SNP_calling/03.gatk/gatk-4.2.1.0/gatk
CPU=36 #12
SCRATCHDIR=`pwd`/tmp
sPath=./gvcf
vcfs="$sPath/MOPhn-1*.gz $sPath/MOPhn-3*.gz $sPath/MOPhn-4*.gz $sPath/MOPhn-5*.gz $sPath/MOPhn-7*.gz $sPath/MOPsg-1*.gz $sPath/MOPsg-3*.gz $sPath/MOPsg-4*.gz $sPath/MOPsg-5*.gz $sPath/MOPsg-6*.gz $sPath/MOPsg-7*.gz $sPath/MHKmls-12*.gz $sPath/MHKmls-13*.gz $sPath/MHKmls-21*.gz $sPath/MHKmls-5*.gz $sPath/MHKmls-7*.gz"
ulimit -n 99999 #change file handler limit!

mkdir -p $SCRATCHDIR



#STEP 6: combine gvcf
STEP=06.combinegvcf.$outstem

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    ls $vcfs > gvcfs.list
    $gatk --java-options "-Xmx180G -XX:ParallelGCThreads=4" CombineGVCFs -R $genome --variant gvcfs.list -O $outstem.combined_gvcf.vcf > $STEP.log 2>&1 \
    && touch  $STEP.done

    bgzip -@ $CPU $outstem.combined_gvcf.vcf >> $STEP.log 2>&1 && tabix -p vcf $outstem.combined_gvcf.vcf.gz >> $STEP.log 2>&1

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi

#STEP 7: genotype gvcf
STEP=07.genotypegvcf.$outstem

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    $gatk --java-options '-Xmx120g' GenotypeGVCFs -O $outstem.genotyped.g.vcf -R $genome -V $outstem.combined_gvcf.vcf.gz  -all-sites true > $STEP.log 2>&1 \
    && touch  $STEP.done

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi

#STEP 8: compress
STEP=08.compress.$outstem

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    bgzip -@ $CPU $outstem.genotyped.g.vcf > $STEP.log 2>&1 && tabix -p vcf $outstem.genotyped.g.vcf.gz >> $STEP.log 2>&1 \
    && touch  $STEP.done

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi
