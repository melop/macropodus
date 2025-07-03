Pop=MOPhn

snpeffjar=/data/projects/zwang/snpEff/snpEff.jar
REF=Macropodus_hongkongensis_2.0

java -Xmx26g -jar $snpeffjar $REF ../joinedcalled.snps.indels.vcf.gz | bgzip -c > $Pop.merged.snpeff.vcf.gz
