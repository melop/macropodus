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
export QUARRY_PATH=/data/software/CodingQuarry_v2.0/QuarryFiles/
export PATH=$PATH:/data/software/GlimmerHMM/bin:/data/software/GlimmerHMM/train/

funannotate iprscan -i trained -m docker -c $CPU
