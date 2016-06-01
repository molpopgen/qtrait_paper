#!sh

#0.1 deleted from both loops right now..
for H2 in 1.0
do
    for OPT in 0.1 0.5 1
    do
	for mu in 0.00025 0.005 0.001
	do
	    qsub hpc/samples_time_additive.sh $mu $H2 H2_"$H2"_OPT_"$OPT"_mu_"$mu".samples.pickle.gz $OPT 10
	done
    done
done
