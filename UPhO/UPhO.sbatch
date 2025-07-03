#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 80
#SBATCH -p blade,himem,hugemem


export PATH=/data/software/UPhO/:$PATH #add the path of UPHO executables

PROT=allproteins.fa #File containing all protein sequences formatted by preprocessAA.php

cat allENSproteins.fa  allmacropodusproteins.fa > $PROT

Blast_helper.sh -i $PROT -c 60 > blastout.log 2>&1



