#!/usr/bin/bash
#SBATCH -p hugemem,himem,gpu,long

module load clusterbasics
module load hyphy
php getNonsyn.php > 1_getNonsyn.log 2>&1
