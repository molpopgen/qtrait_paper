#!/bin/bash

for i in mlocus_pickle/*.tar
do
	echo $i
	n=`basename $i`
	qsub hpc/mlocus11_process_pop_var.sh $n
done

