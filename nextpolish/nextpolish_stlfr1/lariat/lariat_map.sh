LARIAT=/data/software/lariat/go/bin/lariat
outdir=lariat_mapped
threads=62

mkdir -p $outdir

$LARIAT -output $outdir  -threads $threads  -genome /data/projects/rcui/mop/nextpolish/3ddna/genome.nextpolish.0.review.fasta  -reads sorted_for_lariat.fq.gz > lariat.log 2>&1

