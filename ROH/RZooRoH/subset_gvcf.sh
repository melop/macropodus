bcftools view -s /data/projects/zwang/m.hk/GATK/mop_on_mhk/DSY-1.bam,/data/projects/zwang/m.hk/GATK/DSY-2.bam,/data/projects/zwang/m.hk/GATK/macro_for_add/DSY-3.bam,/data/projects/zwang/m.hk/GATK/DSY-4.bam,/data/projects/zwang/m.hk/GATK/DSY-5.bam,/data/projects/zwang/m.hk/GATK/macro_for_add/DSY-6.bam,/data/projects/zwang/m.hk/GATK/DSY-7.bam,mhk,/data/projects/zwang/m.hk/GATK/SZ-2.bam,/data/projects/zwang/m.hk/GATK/SZ-3.bam,/data/projects/zwang/m.hk/GATK/SZ-4.bam,MHKqns-1,MHKqns-2,MHKqns-3  -Oz -o MHKmlh_qns_MOPdsy.gvcf.gz allpops.g.vcf.gz
bcftools index MHKmlh_qns_MOPdsy.gvcf.gz

