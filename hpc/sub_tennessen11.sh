#!/bin/bash

for mu in 0.00025 # 0.001 0.005
do
	for opt in 1.0 # 0 0.1 0.5 1.0
	do
		qsub hpc/tennessen11.sh $mu $opt $RANDOM pop.mu"$mu".opt"$opt"
	done
done
