#!/bin/bash

#$ -q krt,krti,krti,bio,sf
#$ -pe openmp 64
#$ -R y

module load krthornt/thorntonlab

cd $SGE_O_WORKDIR
sampler=$1
mutrate=$2
ofile=$3
tsample=$4
optimum=$5
VS=$6
sigmu=$7
seed=$8
echo "seed = $seed"
#This will default to 10 loci, 64 cores, 16 batches
#Args 9 and 10 should be --fixations filename to record fixation times & other data
#Arg 11 will be the sample size for case when sampler == popgen
if [ -z ${11+x} ]
then
    #There is no arg 11
    /usr/bin/time -f "%e %M" -o multi_region.$sampler.$optimum.$VS.$mutrate.time python python/multi_region.py --sampler $sampler -m $mutrate -F $ofile -t $tsample -O $optimum -S $VS -e $sigmu -N 10000 -s $seed $9 ${10}
else
    echo "python python/multi_region.py --sampler $sampler -m $mutrate -F $ofile -t $tsample -O $optimum -S $VS -e $sigmu -N 10000 -s $seed $9 ${10} --nsam ${11}"
    /usr/bin/time -f "%e %M" -o multi_region.$sampler.$optimum.$VS.$mutrate.time python python/multi_region.py --sampler $sampler -m $mutrate -F $ofile -t 1000 --t2 50 -O $optimum -S $VS -e $sigmu -N 10000 -s $seed $9 ${10} --nsam ${11}
fi
