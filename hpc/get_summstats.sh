#!/bin/bash

#$ -q krti
#$ -t 1-36

cd $SGE_O_WORKDIR

module load krthornt/thorntonlab

INFILE=`ls *.samples.pickle.gz|head -n $SGE_TASK_ID|tail -n 1`
OFILE=`basename $INFILE .pickle.gz`
python python/get_summstats.py $INFILE 16 $OFILE.stats.h5