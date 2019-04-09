#!/bin/bash

#$ -q krt2
#$ -pe openmp 64-128

cd $SGE_O_WORKDIR
module load krthornt/anaconda
source activate fwdpy11_0_3_0

dir=$1
sdir=$2
SEED=$3

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 analyze_ancient_samples.py $dir/$sdir 1 256 \
$dir"_"$sdir"_gscan.sqlite3" $dir"_"$sdir"_qtrait.sqlite3" $dir"_"$sdir"_ld.sqlite3" $SEED $CORES
