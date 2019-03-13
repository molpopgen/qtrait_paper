#!/bin/bash

qsub -N MAIN_LOMU run.sh 0.00025 0.25 lomu/replicate lomu_seeds.txt
qsub -N MAIN_MIDMU run.sh 1e-3 0.25 midmu/replicate midmu_seeds.txt
qsub -N MAIN_HIMU run.sh 5e-3 0.25 himu/replicate himu_seeds.txt
