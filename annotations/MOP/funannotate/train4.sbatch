#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 64
#SBATCH --mem=100G
source /data/software/funannotate/env.sh

GENOME=/data/projects/rcui/mop/annotations/repeatmasker/scf.softmasked.fa
TRINITY_TRANSCRIPT=all.transcripts.fa
SPECIES="Macropodus opercularis"
STRAIN="DSY2021" #di shui yan park
CPU=64

funannotate train -i $GENOME -o trained  \
--left /data/projects/rcui/mop/annotations/funannotate/train/trained/training.bak/left.fq.gz \
--right /data/projects/rcui/mop/annotations/funannotate/train/trained/training.bak/right.fq.gz  \
--no_trimmomatic \
--no_normalize_reads \
--pasa_min_avg_per_id 85 \
--jaccard_clip \
--memory 100G \
--trinity "$TRINITY_TRANSCRIPT" \
--species "$SPECIES" \
--pasa_db mysql \
--max_intronlen 200000 \
--aligners minimap2 minimap2 \
--strain "$STRAIN" --cpus $CPU 
