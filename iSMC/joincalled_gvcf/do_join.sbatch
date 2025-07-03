#!/bin/bash
#SBATCH -p long,blade,himem
#SBATCH -c 36
#SBATCH --mem 120G

#do gatk.sh for each sample.bam file
module load bwa
module load java
module load samtools

bash joined_call.sh
