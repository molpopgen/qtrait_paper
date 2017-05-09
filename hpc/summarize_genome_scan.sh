#!/bin/bash

#$ -q krt,krti

opt=$1
mu=$2

module load krthornt/anaconda

cd $SGE_O_WORKDIR

python ../../python/get_summstats.py $opt $mu
