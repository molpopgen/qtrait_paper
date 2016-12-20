#!/bin/bash

#$ -pe openmp 64
#$ -q krt,krti,pub64,free64,free72i

cd $SGE_O_WORKDIR

module load krthornt/thorntonlab

python python/validate_additive.py
