#!/bin/bash

for timing in pre post
do
	for sel in neutral
	do
		NFILES=`find genome_scan -maxdepth 1 -name "*$sel.$timing"_shift"*.gz" | wc -l`
		SCRIPTNAME=genome_scan_summstats.$timing.$sel.sh
		if true ; then
		cat <<- EOF > $SCRIPTNAME
		#!/bin/bash
		#\$ -q krt,krti,bio,pub64,sf
		#\$ -pe openmp 8
		#\$ -t 1-$NFILES
		cd \$SGE_O_WORKDIR

		module load krthornt/thorntonlab

		infile=\`find genome_scan -maxdepth 1 -name "*$sel.$timing"_shift"*.gz" | head -n \$SGE_TASK_ID | tail -n 1\`
		stub=\`basename \$infile .gz\`
		outfile=genome_scan/summstats/\$stub.summstats.gz

		./src/summstats \$infile \$outfile 0.05 0.1 \$CORES
		EOF
		qsub -N stats.$timing.$sel $SCRIPTNAME
	fi
	done
done
