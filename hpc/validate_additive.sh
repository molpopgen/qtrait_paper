#!/bin/bash

#$ -pe openmp 64
#$ -q krt

cd $SGE_O_WORKDIR

module load krthornt/thorntonlab

python python/validate_additive.py