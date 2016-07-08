#!/bin/bash

VS=1
mu=1e-3
sigmu=0.25

qsub -N HIMU2 hpc/run_single_region.sh stats $mu himu2.h5 100 0.5 $VS $sigmu $RANDOM

#Drop mutation rate 10x, up VS to keep 4muVS 
#but do NOT hold mutational variance constant
VS=10
mu=1e-4

qsub -N LOMU2 hpc/run_single_region.sh stats $mu lomu2.h5 100 0.5 $VS $sigmu $RANDOM

#Do it again...
VS=100
mu=1e-5

qsub -N VLOMU2 hpc/run_single_region.sh stats $mu vlomu2.h5 100 0.5 $VS $sigmu $RANDOM
