#!/bin/bash

for mu in 0.005
do
	for opt in 1.0 0.5 0.1
	do
		qsub -q krti,krt,bio hpc/mlocus_plus_neutral11.sh $mu $opt $RANDOM pop.mu"$mu".opt"$opt"
	done
done
for mu in 0.001 0.00025 
do
	for opt in 1.0 0.5 0.1
	do
		qsub hpc/mlocus_plus_neutral11.sh $mu $opt $RANDOM pop.mu"$mu".opt"$opt"
	done
done
