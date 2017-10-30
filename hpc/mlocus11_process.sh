#!/bin/bash

#$ -q krt,krti,bio
#$ -pe openmp 16-32

cd $SGE_O_WORKDIR/mlocus_pickle

module load krthornt/anaconda/3

INFILE=$1
/usr/bin/time -f "%e %M" -o time.out python ../python/mlocus11_process.py -t $INFILE --nprocs $CORES
