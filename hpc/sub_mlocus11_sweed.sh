#!/bin/bash

NJOBS=1
for i in mlocus_pickle/pop*.tar
do
    n=`basename $i .tar`
    if [ ! -e mlocus_pickle/$n"_sweed.db" ]
    then
        echo $i $n
        NAME="SWEEDJOB_"$NJOBS
        if [ $NJOBS -lt 4 ]
        then
            echo $NAME
            qsub -N $NAME hpc/run_mlocus11_sweed.sh $n.tar $n"_sweed.db"
        else
            X=$(($NJOBS-3))
            LASTJOB="SWEEDJOB_"$X
            echo $NAME $LASTJOB
            qsub -N $NAME -hold_jid $LASTJOB hpc/run_mlocus11_sweed.sh $n.tar $n"_sweed.db"
        fi
    fi
    NJOBS=$(($NJOBS+1))
done
