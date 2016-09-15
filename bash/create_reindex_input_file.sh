#!/bin/bash

for i in *_popgen.h5
do
	n=`basename $i .h5`
	echo $i reindexed/$n.reindexed.h5 summstats
done

for i in *_stats.h5
do
	n=`basename $i .h5`
	echo $i reindexed/$n.reindexed.h5 stats
done

