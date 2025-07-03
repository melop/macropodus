#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 64
#SBATCH --mem=100G
source /data/software/funannotate/env.sh

GENOME=/data/projects/rcui/mop/annotations/repeatmasker/scf.softmasked.fa
STRINGTIE=/data/projects/rcui/mop/annotations/funannotate/rna_evidence/stringtie/longreads_assm.out.gtf
SPECIES="Macropodus opercularis"
STRAIN="DSY2021" #di shui yan park
CPU=64
export QUARRY_PATH=/data/software/CodingQuarry_v2.0/QuarryFiles/
export PATH=$PATH:/data/software/GlimmerHMM/bin:/data/software/GlimmerHMM/train/

 funannotate predict -i $GENOME \
            -o trained -s "$SPECIES" --strain $STRAIN --cpus $CPU \
--max_intronlen 200000 \
--weights augustus:10 hiq:10 pasa:10 snap:0 genemark:0 glimmerhmm:0 proteins:10 transcripts:1 \
--optimize_augustus \
--busco_db actinopterygii \
--busco_seed_species zebrafish \
--organism other \
--protein_evidence /data/projects/rcui/mop/annotations/funannotate/prot_evidence/all.proteins.fa \
--stringtie $STRINGTIE \
--keep_no_stops \
--min_protlen 10 \
--repeat_filter none
#--repeats2evm \
#--no-evm-partitions \
# StartWeights = {'augustus': 1, 'hiq': 2, 'genemark': 1, 'pasa': 6,'codingquarry': 0, 'snap': 1, 'glimmerhmm': 1, 'proteins': 1, 'transcripts': 1}
