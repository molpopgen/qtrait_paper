#!/bin/bash

#$ -q krt2
#$ -pe openmp 128

cd $SGE_O_WORKDIR/mlocus_pickle

module load krthornt/anaconda/3

INFILE=$1
python ../python/mlocus11_summarize_population_variation.py -t $INFILE --nprocs $CORES


