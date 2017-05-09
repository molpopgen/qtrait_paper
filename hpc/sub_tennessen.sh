#!/bin/bash

#for mu in 0.00025 0.001 0.005
for mu in 0.005
do
	#for opt in 0 0.1 0.5 1.0
	for opt in 1.0	
	do
		qsub hpc/tennessen.sh $mu $opt $RANDOM stats.mu"$mu".opt"$opt".h5
	done
done
