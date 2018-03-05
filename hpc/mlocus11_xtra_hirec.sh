#!/bin/bash

#$ -q krt,krti,bio
#$ -pe openmp 64-72

module load krthornt/anaconda/3

cd $SGE_O_WORKDIR

mu=$1
opt=$2
seed=$3
plarge=$4
vsopt=$5
db1=$6
db2=$7

python python/multi_region11_explore.py --ncores $CORES --nreps 256 --mu $mu --opt $opt --seed $seed --plarge $plarge --vsopt \
    $vsopt --statsdb $db1 --scandb $db2 --rho 1000.0

