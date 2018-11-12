#!/bin/bash

for mu in 0.005 0.001 0.00025
do
    for plarge in 0.05 0.1 0.25 0.5 0.75
    do
        SUFFIX=mu$mu".plarge"$plarge".opt1.sqlite3"
        echo $SUFFIX
        qsub run_gamma.sh $plarge $mu qtrait_gamma.$SUFFIX gscan_gamma.$SUFFIX fixations_gamma.$SUFFIX
    done
done
