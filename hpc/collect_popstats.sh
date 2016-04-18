#!/bin/bash

#$ -q krt,krti
#Not sure how much RAM this needs.
# $ -pe openmp 64  

cd $SGE_O_WORKDIR

module load krthornt/thorntonlab

/usr/bin/time -f "%e %M" -o collect_popstats.time.txt python python/collect_popstats.py