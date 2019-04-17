#!/bin/bash

I=1
# for dir in main_results lowrec # hirec
for dir in hirec
do
    for sdir in lomu midmu himu
    do
        SEED=`head -n $I analyze_seeds.txt | tail -n 1`
        qsub analyze.sh $dir $sdir $SEED
        I=$(($I+1))
    done
done
