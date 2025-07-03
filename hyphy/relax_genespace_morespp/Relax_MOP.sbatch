#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/RelaxNothos_slurm-%a.out
#SBATCH -e slurmlog/RelaxNothos_slurm-%a.err
#SBATCH -a 0-199
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

#SET THIS TO THE CORRECT NUMBER
TreeFolder="../exportAln_genespace_morespp/Relax_MOP"
OutDIR=`basename $TreeFolder`;

TOTALPARTS=100

mkdir -p $OutDIR
mkdir -p logs

for i in `seq 0 99`; do
php hyphyrelax.php -N ${TOTALPARTS} -f ${i} -T $TreeFolder -o $OutDIR > logs/$i.log 2>&1 &

done
wait
exit

