#!/bin/bash

#$ -q krt,krti,bio
#$ -pe openmp 64-72

module load krthornt/anaconda/3

cd $SGE_O_WORKDIR

mu=$1
opt=$2
seed=$3
picklefile=$4

/usr/bin/time -f "%e %M" -o mlocus_plus_neutral_pickle/time.txt python python/multi_region_plus_neutral11.py --ncores $CORES --nreps 1024 --mu $mu --opt $opt --seed $seed --stub mlocus_plus_neutral_pickle/$picklefile --tarfile mlocus_plus_neutral_pickle/pops.mu$mu.opt$opt.tar
