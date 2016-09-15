#!/bin/bash

for i in reindexed/*.h5
do
	n=`basename $i .reindexed.h5`
	if [ -e $n.h5 ]
	then
		mv -f $i $n.h5
	fi
done
