#!/bin/bash

#$ -q krt,krti,bio,sf
#$ -pe openmp 64-72
#$ -ckpt blcr

module load krthornt/anaconda/3

cd $SGE_O_WORKDIR

mu=$1
opt=$2
seed=$3
picklefile=$4

/usr/bin/time -f "%e %M" -o mlocus_pickle/time.txt python python/multi_region11.py --ncores $CORES --nreps 256 --mu $mu --opt $opt --seed $seed --stub mlocus_pickle/$picklefile --tarfile mlocus_pickle/pops.mu$mu.opt$opt.tar
