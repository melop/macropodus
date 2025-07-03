
mkfifo all.fulllength.ont.cdna.fq
cat ../../isoseq/*.full_length*.fq >  all.fulllength.ont.cdna.fq &
REF=/data/projects/rcui/mhk/nextdenovo/nextpolish_stlfr2/mhk.LG.nextdenovo.sgs.2.sorted.fa
CPU=64

 minimap2 -ax splice -t $CPU --cs -u b -G 1000000 $REF  all.fulllength.ont.cdna.fq | samtools sort --reference $REF -@ 12 -o mhk.ont.cdna.alignments.bam - > aln.log 2>&1
stringtie -p $CPU -L -o longreads_assm.out.gtf  mhk.ont.cdna.alignments.bam > stringtie.log 2>&1
/data/software/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl longreads_assm.out.gtf $REF > mhk.assembled_transcripts.fa
