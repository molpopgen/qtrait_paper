#!/bin/bash

#$ -q krt,krti,bio,pub64
#$ -pe openmp 64
#$ -R y

module load krthornt/thorntonlab

cd $SGE_O_WORKDIR
SEED=$6
python python/ten_loci_samples_additive.py -m $1 -H $2 -F $3 -O $4 -s $SEED -t $5
