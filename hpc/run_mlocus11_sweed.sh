#!/bin/bash

#$ -q krt,krti,bio
#$ -pe openmp 64-72

module load krthornt/anaconda

cd $SGE_O_WORKDIR/mlocus_pickle

INFILE=$1
OUTFILE=$2
MU=`echo $INFILE | sed 's/^.*u//' | sed 's/\.opt.*//'`
OPT=`echo $INFILE | sed 's/^.*opt//' | sed 's/\.tar//'`
echo "OPT=$OPT"
P=`echo "$CORES/2" | bc`
python ../python/mlocus11_sweed.py --tarfile $INFILE --outfile $OUTFILE --mu $MU --opt $OPT --nprocs $P
