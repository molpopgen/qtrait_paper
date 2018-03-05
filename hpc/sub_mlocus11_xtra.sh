#!/bin/bash

for mu in 0.001 # 0.00025 0.005
do
	for opt in 1.0
    do
        for plarge in 0.0 0.01 0.10
        do
            for vsopt in 1 0.5 0.1
            do
                db1=mu"$mu".opt"$opt".plarge"$plarge".vsopt"$vsopt".stats.sqlite3
                db2=mu"$mu".opt"$opt".plarge"$plarge".vsopt"$vsopt".scan.sqlite3
                qsub -q krti,krt,bio hpc/mlocus11_xtra.sh $mu $opt $RANDOM $plarge $vsopt mlocus11_xtra_sims/$db1 mlocus11_xtra_sims/$db2
            done
        done
    done
done
