#!/bin/bash

for i in mlocus_pickle/*.tar
do
	echo $i
	n=`basename $i`
	qsub hpc/mlocus11_gamma.sh $n
done
