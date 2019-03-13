#!/bin/bash

qsub -N LOREC_LOMU run.sh 0.00025 0.25 lomu/replicate lomu_seeds.txt
qsub -N LOREC_MIMDU run.sh 1e-3 0.25 midmu/replicate midmu_seeds.txt
qsub -N LOREC_HIMU run.sh 5e-3 0.25 himu/replicate himu_seeds.txt
