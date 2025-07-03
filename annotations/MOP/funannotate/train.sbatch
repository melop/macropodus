#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 64
#SBATCH --mem=100G
source /data/software/funannotate/env.sh

GENOME=/data/projects/rcui/mop/annotations/repeatmasker/scf.softmasked.fa
TRINITY_TRANSCRIPT=/data/projects/rcui/mop/trinity/transcoder/Trinity.cdhitout.filtered.fa
NANOPORE=`realpath /data/projects/rcui/mhk/annotations/stringtie_isoseq/mhk.assembled_transcripts.fa`
SPECIES="Macropodus opercularis"
STRAIN="DSY2021" #di shui yan park
CPU=64

funannotate train -i $GENOME -o trained  \
--left /data/projects/shareddata/mop.stLFR/RNA/NGS/trim/brain.paired_1.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/fin.paired_1.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/gill.paired_1.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/intestine.paired_1.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/liver.paired_1.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/muscle.paired_1.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/skin.paired_1.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/snout.paired_1.fq.gz \
--right /data/projects/shareddata/mop.stLFR/RNA/NGS/trim/brain.paired_2.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/fin.paired_2.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/gill.paired_2.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/intestine.paired_2.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/liver.paired_2.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/muscle.paired_2.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/skin.paired_2.fq.gz \
/data/projects/shareddata/mop.stLFR/RNA/NGS/trim/snout.paired_2.fq.gz \
--no_trimmomatic \
--memory 100G \
--trinity "$TRINITY_TRANSCRIPT" \
--nanopore_cdna "$NANOPORE" \
--species "$SPECIES" \
--pasa_db mysql \
--max_intronlen 1000000 \
--strain "$STRAIN" --cpus $CPU

