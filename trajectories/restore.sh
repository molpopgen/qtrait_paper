#!/bin/bash

#$ -q krti,krt
#$ -pe openmp 32-72
#$ -t 1-9

cd $SGE_O_WORKDIR
module load krthornt/anaconda/3

INFILE=`ls *_dump.gz | head -n $SGE_TASK_ID | tail -n 1`
OFILE=`basename $INFILE _dump.gz`

if [ ! -e $OFILE.db ]
then
	zcat $INFILE | sqlite3 $OFILE.db
fi
