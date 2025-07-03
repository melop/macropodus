#!/usr/bin/bash
#SBATCH -c 40
#SBATCH -p long
#SBATCH --mem=100G

ASSEMBLY=/data/projects/rcui/mop/nextpolish/nextpolish_stlfr2/nextpolish/mop.LG.nextdenovo.sgs.2.fa # assembly to mask
REPBASELIB=`realpath ./all_te.fa` #repeat library, with killifish repeats
THREADS=40

ln -sf $ASSEMBLY scf.fa


singularity run docker://dfam/tetools:latest RepeatMasker -gff -html -pa $THREADS -lib $REPBASELIB  scf.fa

