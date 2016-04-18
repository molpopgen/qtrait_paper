#!/bin/bash

#$ -q krt,krti,krti,bio,pub64,sf
#$ -pe openmp 64
#$ -R y

module load krthornt/thorntonlab

cd $SGE_O_WORKDIR
python python/trajectories_time_additive.py -m $1 -H $2 -O $3 -s $4 --fixed $5 --ages $6 --traj $7
