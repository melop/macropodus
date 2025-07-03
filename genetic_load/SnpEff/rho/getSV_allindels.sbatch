#!/usr/bin/bash
#SBATCH -p hugemem,himem,gpu,long

module load clusterbasics
module load hyphy
php getSV_allindels.php > 2_getSV_allindels.log 2>&1
