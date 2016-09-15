#!/bin/bash

#$ -pe openmp 2
#$ -q krt,bio,krti,pub64
#$ -t 1-27

module load krthornt/thorntonlab

cd $SGE_O_WORKDIR

ifile=`head -n $SGE_TASK_ID h5list.txt|tail -n 1 | cut -d" " -f1`
ofile=`head -n $SGE_TASK_ID h5list.txt|tail -n 1 | cut -d" " -f2`
stat=`head -n $SGE_TASK_ID h5list.txt|tail -n 1 | cut -d" " -f3`

python python/reindexh5.py $ifile $ofile $stat

