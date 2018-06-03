#!/bin/bash

#$ -q krt,krti,bio
#$ -pe openmp 32

cd $SGE_O_WORKDIR/mlocus_pickle

module load krthornt/anaconda/3

INFILE=$1
python ../python/mlocus11_omega_max.py -t $INFILE --nprocs 8


