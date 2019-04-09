#!/bin/bash

#$ -q krt2
#$ -pe openmp 32-128

cd $SGE_O_WORKDIR
module load krthornt/anaconda
source activate fwdpy11_0_3_0

python3 create_fixations_database.py $1 1 256 $2 $CORES
