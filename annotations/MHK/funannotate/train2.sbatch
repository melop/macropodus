#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 64
#SBATCH --mem=100G
source /data/software/funannotate/env.sh

GENOME=/data/projects/rcui/mhk/annotations/repeatmasker/scf.softmasked.fa
TRINITY_TRANSCRIPT=/data/projects/rcui/mop/annotations/funannotate/rna_evidence/all.transcripts.fasta #/data/projects/rcui/mop/trinity/transcoder/Trinity.cdhitout.filtered.fa
#NANOPORE=`realpath /data/projects/rcui/mhk/annotations/stringtie_isoseq/mhk.assembled_transcripts.fa`
SPECIES="Macropodus hongkongensis"
STRAIN="SZDP2020" #Shenzhen Daping bay
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
--max_intronlen 1000000 \
--aligners minimap2 gmap \
--strain "$STRAIN" --cpus $CPU

#--nanopore_cdna "$NANOPORE" 
