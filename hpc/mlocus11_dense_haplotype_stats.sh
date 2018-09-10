#!/bin/bash

#$ -q bsg2,krti
#$ -pe openmp 64-128
#$ -ckpt blcr

module load krthornt/anaconda/3
source activate fwdpy11_0_2_0

cd $SGE_O_WORKDIR

mu=$1
opt=$2
seed=$3
statfile=$4
fixationsfile=$5

python3 python/multi_region11_dense_sampling_haplotype_stats.py --ncores $CORES --nreps 256 --mu $mu \
--opt $opt --seed $seed --outfile $statfile --fixationsfile $fixationsfile --simlen 15 \
--nsam 50

