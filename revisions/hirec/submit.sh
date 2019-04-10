#!/bin/bash

#qsub -N HIREC_LOMU run.sh 0.00025 0.25 lomu/replicate lomu_seeds.txt
#qsub -N HIREC_MIMDU run.sh 1e-3 0.25 midmu/replicate midmu_seeds.txt
qsub -N HIREC_HIMU run.sh 5e-3 0.25 himu/replicate himu_seeds.txt
