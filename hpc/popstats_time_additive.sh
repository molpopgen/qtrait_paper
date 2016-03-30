#!/bin/bash

#$ -q krt,krti,bio,pub64,sf
#$ -pe openmp 64
#$ -R y

module load krthornt/thorntonlab

cd $SGE_O_WORKDIR
SEED=$RANDOM
python python/popstats_time_additive.py -m $1 -H $2 -F $3 -O $4 -s $SEED -t $5
