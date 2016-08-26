#!sh
ISEED=202
I=1
nsam=500
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
	    ofile=H2_"$H2"_OPT_"$OPT"_mu_"$mu"_sigmu0.25_10regions_trajectories.h5
	    #Skip sf queue b/c of possible RAM issues...
	    qsub -N MTRAJ -q krt,krti,krti,bio hpc/run_multi_region.sh freq $mu $ofile 1 $OPT 1 0.25 $SEED 
	    #qsub hpc/mloc_samples_time_additive.sh $mu $H2 H2_"$H2"_OPT_"$OPT"_mu_"$mu".mloc_samples.h5 $OPT 10 $SEED
	done
    done
done
