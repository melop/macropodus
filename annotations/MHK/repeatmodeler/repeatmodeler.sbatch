#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 48
#SBATCH --mem=250G

ASSEMBLY=/data/projects/rcui/mhk/nextdenovo/nextpolish_stlfr2/mhk.LG.nextdenovo.sgs.2.fa # the assembly to look at
SP=MHK
REPMODEL="singularity run docker://dfam/tetools:latest " #command of repeatmodeler
THREADS=48

ln -sf $ASSEMBLY scf.fa

$REPMODEL BuildDatabase -name $SP -engine ncbi scf.fa

$REPMODEL RepeatModeler -LTRStruct -engine ncbi -pa $THREADS -database $SP
