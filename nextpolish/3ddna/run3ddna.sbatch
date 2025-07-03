#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 40

Contigs=/data/projects/rcui/mop/nextpolish/01_rundir/genome.nextpolish.fasta
MNF=/data/projects/rcui/mop/nextpolish/juicer/aligned/merged_nodups.txt
MAPQ=30
GAPSIZE=1000


/data/software/3d-dna-201008/run-asm-pipeline.sh -m haploid -i 100 -r 3 -q $MAPQ --sort-output -g $GAPSIZE $Contigs $MNF
