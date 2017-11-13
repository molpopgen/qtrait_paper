#!/bin/bash

#$ -q krti
#$ -pe openmp 72

module load krthornt/anaconda

cd $SGE_O_WORKDIR

python mlocus11_sweed_nulldist.py
