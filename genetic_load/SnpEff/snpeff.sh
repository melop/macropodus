Pop=allpops

snpeffjar=/data/projects/zwang/snpEff/snpEff.jar
REF=Macropodus_hongkongensis_2.0

java -Xmx48g -jar $snpeffjar $REF ../allpops.withoutgroup.norm.vcf.gz | bgzip -c > $Pop.merged.snpeff.vcf.gz
