nextpolish="python /data/software/NextPolish/lib/nextpolish1.py"
threads=24

lariatbam=`realpath ..`/lariat/lariat_mapped/bc_sorted_bam.bam
input=/data/projects/rcui/mop/nextpolish/3ddna/genome.nextpolish.0.review.fasta 

samtools view --threads 12 -F 0x4 -u $lariatbam  |samtools fixmate -m --threads 12  - -| samtools sort -m 2g --threads 24 -|samtools markdup --threads 24 -r - sgs.sort.bam

samtools index -@ ${threads} sgs.sort.bam;

$nextpolish -g ${input} -t 1 -p ${threads} -s sgs.sort.bam > genome.nextpolish.stlfr.1.fa;
funannotate sort -i genome.nextpolish.stlfr.1.fa -o mop.LG.nextdenovo.sgs.1.fa -b mopscf

