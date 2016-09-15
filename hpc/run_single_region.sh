#!/bin/bash

#$ -q krt,krti,krti,bio,pub64,sf
#$ -pe openmp 64
#$ -R y

module load krthornt/thorntonlab

cd $SGE_O_WORKDIR
sampler=$1
mutrate=$2
ofile=$3
tsample=$4
optimum=$5
VS=$6
sigmu=$7
seed=$8
python python/single_region.py --sampler $sampler -m $mutrate -F $ofile -t $tsample -O $optimum -S $VS -e $sigmu -N 5000 -s $seed
