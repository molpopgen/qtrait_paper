#!/bin/bash

#$ -q krt,krti,bio
#$ -pe openmp 8-16

cd $SGE_O_WORKDIR/mlocus_pickle

module load krthornt/anaconda/3

INFILE=$1
mu=`echo $INFILE | cut -d "." -f2,3| sed 's/[a-z]//g'`
opt=`echo $INFILE | cut -d "." -f4,5| sed 's/[a-z]//g'`
ofile=`basename $INFILE .tar`
echo "$INFILE $mu $opt $ofile $CORES"
/usr/bin/time -f "%e %M" -o time.out python ../python/mlocus11_gamma_hat.py -t $INFILE --nprocs $CORES --mu $mu --opt $opt --outfile $ofile.gamma_threshold.gz
