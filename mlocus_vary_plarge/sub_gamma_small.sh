#!/bin/bash

for mu in 0.005 0.001 0.00025
do
    for plarge in 0.05 0.5 0.9
    do
        SUFFIX=mu$mu".plarge"$plarge".opt1.sqlite3"
        echo $SUFFIX
        qsub run_gamma_small.sh $plarge $mu qtrait_gamma_small.$SUFFIX gscan_gamma_small.$SUFFIX fixations_gamma_small.$SUFFIX
    done
done

