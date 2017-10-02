#!/bin/bash

#$ -q krt,krti,bio
#$ -pe openmp 8-16

cd $SGE_O_WORKDIR/mlocus_pickle

module load krthornt/anaconda/3

INFILE=$1
python ../python/mlocus11_fixations.py -t $INFILE --nprocs $CORES
