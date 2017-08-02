#!/bin/bash

#$ -q krti
#$ -pe openmp 32

cd $SGE_O_WORKDIR

module load krthornt/anaconda

Rscript process_frequency_data.R *_merged.db
