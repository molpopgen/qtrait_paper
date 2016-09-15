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
	    if [ $I -eq 3 ]
		then
		    echo $SEED
		    echo $I
		    ofile=H2_"$H2"_OPT_"$OPT"_mu_"$mu"_sigmu0.25_traj
		    qsub -N FATES -q krt,krti,bio hpc/run_single_region.sh freq $mu $ofile 1 $OPT 1 0.25 $SEED
	    	#qsub hpc/fates_additive.sh $mu $H2 $OPT $SEED  H2_"$H2"_OPT_"$OPT"_mu_"$mu".ages.h5 H2_"$H2"_OPT_"$OPT"_mu_"$mu".traj.h5
	    fi	

		    I=$(($I+1))
	done
    done
done
