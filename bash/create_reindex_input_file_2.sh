
#!/bin/bash

for i in *_popgen_fixations.h5
do
	n=`basename $i .h5`
	echo $i reindexed/$n.reindexed.h5 fixations
done
