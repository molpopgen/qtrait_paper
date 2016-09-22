#!/bin/bash

#SweepFinder is annoying.
#It requires an ms-like header at the start of input.
#
#It also makes a tempfile w/o using the proper conventions.
#The tempfile means you cannot run multiple processes in same working 
#directory.  If you do that, you get either nonsensical output or a crash 
#with a core dump.
#
#It also cannot deal with fixed sites in the data.
#To deal with this:
#
#1. Use command grouping and a subprocess shell to create the ms header.
#2. Use trimFixations to read the original data, removed invariant positions, and write new data 
#to stdout.
#3. Do all the work in a temporary directory for each file
#
for timing in post
do
	for sel in neutral
	do
		NFILES=`find genome_scan -maxdepth 1 -name "*$sel.$timing"_shift"*.gz" | wc -l`
		SCRIPTNAME=genome_scan_sweepfinder.$timing.$sel.sh
		if true ; then
		cat <<- EOF > $SCRIPTNAME
		#!/bin/bash
		#\$ -q krt,krti,bio,pub64,sf
		#\$ -t 1-$NFILES
		cd \$SGE_O_WORKDIR

		module load krthornt/thorntonlab
		module load krthornt/SweepFinder/dev

		infile=\`find genome_scan -maxdepth 1 -name "*$sel.$timing"_shift"*.gz" |sort| head -n \$SGE_TASK_ID | tail -n 1\`
		echo \$infile
		stub=\`basename \$infile .gz\`
		outfile=\$stub.sweepfinder.txt
		mkdir genome_scan/sweepfinder/\$stub
		cd genome_scan/sweepfinder/\$stub
		NREPS=\`zgrep -c segs \$infile\`
		SweepFinder -m 100 <({ echo \"ms 100 \$NREPS -t 100\" & ../../../src/trimFixations ../../\$stub.gz; }) \$outfile
		gzip \$outfile
		mv \$outfile.gz ..
		cd ..
		rm -rf \$stub
		EOF
		qsub -N sf.$timing.$sel $SCRIPTNAME
	fi
	done
done
