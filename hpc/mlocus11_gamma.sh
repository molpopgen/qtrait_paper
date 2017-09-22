#!/bin/bash

#$ -q krt,krti,bio
#$ -pe openmp 64-72

cd $SGE_O_WORKDIR/mlocus_pickle

module load krthornt/anaconda/3

INFILE=$1
mu = `echo $INFILE | cut -d "." -f2,3| sed 's/[a-z]//g'`
opt = `echo $INFILE | cut -d "." -f4,5| sed 's/[a-z]//g'`
/usr/bin/time -f "%e %M" -o time.out python ../python/mlocus11_process.py -t $INFILE --nprocs $CORES --mu $mu --opt $opt
