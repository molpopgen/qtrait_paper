#!/bin/bash

for OPT in 0.1 0.5 1
do
	for mu in 0.00025 0.005 0.001
	do
		qsub ../../hpc/summarize_genome_scan.sh $OPT $mu
	done
done
