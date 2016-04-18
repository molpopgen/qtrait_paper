#!sh

#0.1 deleted from both loops right now..
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
	    qsub hpc/fates_additive.sh $mu $H2 $OPT $SEED H2_"$H2"_OPT_"$OPT"_mu_"$mu".fixations.h5 H2_"$H2"_OPT_"$OPT"_mu_"$mu".ages.h5 H2_"$H2"_OPT_"$OPT"_mu_"$mu".traj.h5
	done
    done
done
