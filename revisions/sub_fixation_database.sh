#!/bin/bash

#for i in main_results lowrec # hirec
for i in hirec
do
    qsub make_fixation_databases.sh $i $i"_fixations.sqlite3"
done

