#!/bin/bash

#$ -q krt2
#$ -pe openmp 64-128

cd $SGE_O_WORKDIR
module load krthornt/anaconda
source activate fwdpy11_0_2_0
PLARGE=$1
MU=$2
QO=$3
GO=$4
FO=$5

python3 ../python/mlocus11_vary_plarge.py --N 5000 --mu $MU \
    --plarge $PLARGE \
    --nsam 50 \
    --nreps 128 \
    --ncores $CORES \
    --opt 1.0 \
    --seed $RANDOM \
    --qoutfile $QO \
    --goutfile $GO \
    --foutfile $FO \
    --gamma 1.0
