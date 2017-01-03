#!/bin/bash

#$ -pe openmp 64
#$ -q bio,krt,krti,pub64
#$ -R y
cd $SGE_O_WORKDIR

module load krthornt/thorntonlab

python python/validate_additive.py
