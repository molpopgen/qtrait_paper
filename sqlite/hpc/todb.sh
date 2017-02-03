#!/bin/bash

#$ -q krt,krti
#$ -pe openmp 8

module load krthornt/anaconda/3
cd $SGE_O_WORKDIR
make -f Makefile.convert -j $CORES
