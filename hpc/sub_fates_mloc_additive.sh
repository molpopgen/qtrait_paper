#!sh

ISEED=202
I=1
for H2 in 1.0
do
    for OPT in 0.1 0.5 1
    do
	for mu in 0.00025 0.005 0.001
	do
	    SEED=`echo "$ISEED*$I"|bc -l`
	    echo $SEED
	    I=$(($I+1))
	    echo $I
	    qsub hpc/mloc_fates_additive.sh $mu $H2 $OPT $SEED H2_"$H2"_OPT_"$OPT"_mu_"$mu".mloc_traj.h5
	done
    done
done
