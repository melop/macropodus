#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 64
#SBATCH --mem=100G
source /data/software/funannotate/env.sh

GENOME=/data/projects/rcui/mhk/annotations/repeatmasker/scf.softmasked.fa
TRINITY_TRANSCRIPT=/data/projects/rcui/mop/trinity/transcoder/Trinity.cdhitout.filtered.fa
NANOPORE=`realpath /data/projects/rcui/mhk/annotations/stringtie_isoseq/mhk.assembled_transcripts.fa`
SPECIES="Macropodus hongkongensis"
STRAIN="SZDP2020" #Shenzhen Daping bay
CPU=64

funannotate train -i $GENOME -o trained  \
--left /data/projects/rcui/mop/annotations/funannotate/train/trained/training/left.fq.gz  \
--right /data/projects/rcui/mop/annotations/funannotate/train/trained/training/right.fq.gz  \
--no_trimmomatic \
--memory 100G \
--trinity "$TRINITY_TRANSCRIPT" \
--nanopore_cdna "$NANOPORE" \
--species "$SPECIES" \
--pasa_db mysql \
--max_intronlen 1000000 \
--strain "$STRAIN" --cpus $CPU

