#!/bin/bash

VS=1
mu=1e-3
sigmu=0.25

qsub -N HIMU hpc/run_single_region.sh stats $mu himu.h5 100 0.5 $VS $sigmu $RANDOM

#Drop mutation rate 10x, up VS to keep 4muVS the same
#increase sigmu to keep mu*simgu^2 the same
VS=10
mu=1e-4
sigmu=0.7905694

qsub -N LOMU hpc/run_single_region.sh stats $mu lomu.h5 100 0.5 $VS $sigmu $RANDOM

#Do it again...
VS=100
mu=1e-5
sigmu=2.5

qsub -N VLOMU hpc/run_single_region.sh stats $mu vlomu.h5 100 0.5 $VS $sigmu $RANDOM
