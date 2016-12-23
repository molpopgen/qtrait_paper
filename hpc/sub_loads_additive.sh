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
	    echo $I
	    ofile=H2_"$H2"_OPT_"$OPT"_mu_"$mu"_sigmu0.25_load.h5
	    qsub -N LOAD -q krt,krti,bio,free64,free72i,pub64,abio -ckpt blcr hpc/run_single_region.sh load $mu $ofile 1 $OPT 1 0.25 $SEED
	    done
	done
done
