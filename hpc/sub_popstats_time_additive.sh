#!sh

ISEED=202
I=1
for H2 in 0.2 #0.1 0.2 0.5
do
    for OPT in 0.1 0.5 1 2.5
    do
	for mu in 0.00025 0.005 0.001
	do
	    SEED=`echo "$ISEED*$I"|bc -l`
	    echo $SEED
	    I=$(($I+1))
	    echo $I
	    qsub hpc/popstats_time_additive.sh $mu $H2 H2_"$H2"_OPT_"$OPT"_mu_"$mu".popstats.h5 $OPT 1 $SEED
	done
    done
done
