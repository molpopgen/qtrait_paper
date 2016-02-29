#!sh

#0.1 deleted from both loops right now..
for H2 in 0.1 0.2 0.5
do
    for OPT in 0.1 0.5 1 5
    do
	for mu in 0.00025 0.005 0.001
	do
	    qsub hpc/fates_additive.sh $mu $H2 $OPT $RANDOM H2_"$H2"_OPT_"$OPT"_mu_"$mu".fixations.h5 H2_"$H2"_OPT_"$OPT"_mu_"$mu".ages.h5
	done
    done
done
