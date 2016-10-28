#!/bin/bash

#$ -q krt,krti,bio
#$ -pe openmp 64-72

module load krthornt/thorntonlab

cd $SGE_O_WORKDIR/tennessen

mu=$1
opt=$2
seed=$3
statfile=$4

echo "Running on $HOSTNAME"
echo "python ../python/TennessenQtraitMloc.py --ncores 32 --nreps 4 --mu $mu --opt $opt --seed $seed --statfile $statfile --tbb $CORES"
python ../python/TennessenQtraitMloc.py --ncores 32 --nreps 4 --mu $mu --opt $opt --seed $seed --statfile $statfile --tbb $CORES
