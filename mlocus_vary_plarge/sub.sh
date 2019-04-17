#!/bin/bash

for mu in 0.005 0.001 0.00025
do
    #for plarge in 0.1 0.5 0.75
    for plarge in 0.75
    do
        SUFFIX=mu$mu".plarge"$plarge".opt1.sqlite3"
        echo $SUFFIX
        qsub run.sh $plarge $mu qtrait.$SUFFIX gscan.$SUFFIX fixations.$SUFFIX
    done
done
