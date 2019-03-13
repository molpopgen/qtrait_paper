#!/bin/bash

#$ -q krt2
#$ -t 1-256
#$ -pe openmp 4

MU=$1
SIGMA=$2
STUB=$3
SEEDFILE=$4

cd $SGE_O_WORKDIR

module load krthornt/anaconda/3
source activate fwdpy11_0_3_0

SEED=`head -n $SGE_TASK_ID $SEEDFILE | tail -n 1`

python ../simulation.py --popsize 5000 --mu $MU --sigma $SIGMA --filename $STUB --seed $SEED  \
    --repid $SGE_TASK_ID --preserve 1 --time 0.04 --num_ind 5000 --rho 10000
