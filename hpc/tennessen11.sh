#!/bin/bash

#$ -q krt,krti,bio
#$ -pe openmp 64-72

module load krthornt/anaconda/3

cd $SGE_O_WORKDIR

mu=$1
opt=$2
seed=$3
picklefile=$4

echo "Running on $HOSTNAME"
echo "python python/TennessenQtraitMloc11.py --ncores 8 --nreps 8 --mu $mu --opt $opt --seed $seed --stub tennessen_pickle/$picklefile"
/usr/bin/time -f "%e %M" -o tennessen_pickle/time.txt python python/TennessenQtraitMloc11.py --ncores 32 --nreps 32 --mu $mu --opt $opt --seed $seed --stub tennessen_pickle/$picklefile --tarfile tennessen_pickle/pops.mu$mu.opt$opt.tar
