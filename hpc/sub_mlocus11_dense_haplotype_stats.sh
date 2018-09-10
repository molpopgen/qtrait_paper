#!/bin/bash

# for mu in 0.005
# do
# 	for opt in 1.0 0.5 0.1
# 	do
# 		qsub -q krti,krt,bio hpc/mlocus11.sh $mu $opt $RANDOM pop.mu"$mu".opt"$opt"
# 	done
# done
# for mu in 0.001 0.00025 
for mu in 0.00025 0.001 0.005
do
    if [ "$mu" == "0.00025" ]
    then
        for opt in 0.5 0.1
        do
            echo "$mu $opt"
            qsub -N dense hpc/mlocus11_dense_haplotype_stats.sh $mu $opt $RANDOM mlocus_pickle/pop.mu"$mu".opt"$opt".dense_hapstats.sqlite3 \
            mlocus_pickle/pop.mu"$mu".opt"$opt".dense_hapstats_fixations.sqlite3
        done
    else
       for opt in 1.0 0.5 0.1
        do
            echo "$mu $opt"
            qsub -N dense hpc/mlocus11_dense_haplotype_stats.sh $mu $opt $RANDOM mlocus_pickle/pop.mu"$mu".opt"$opt".dense_hapstats.sqlite3 \
            mlocus_pickle/pop.mu"$mu".opt"$opt".dense_hapstats_fixations.sqlite3
        done
    fi
done

