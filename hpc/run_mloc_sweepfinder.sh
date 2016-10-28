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
#
#Further, it writes output continuously, and thus is 
#subject to data loss due to NSF lag on a cluster.

#To deal with this:
#
#1. Use command grouping and a subprocess shell to create the ms header.
#2. Use trimFixations to read the original data, removed invariant positions, and write new data 
#to stdout.
#3. Do all the work in a temporary directory for each file
#4. Do all the work on /bio, our BeeGFS partition.  This script must be submitted from /bio/sweepfinder

for timing in post
do
	for sel in neutral
	do
		NFILES=`wc -l /share/kevin2/krthornt/qtrait_paper/genome_scan/neutral_post.txt|cut -f 1 -d" "`
		SCRIPTNAME=genome_scan_sweepfinder.$timing.$sel.sh
		if true ; then
		cat <<- EOF > $SCRIPTNAME
		#!/bin/bash
		#\$ -q krt,krti,bio,pub64,sf,pub8i
		#\$ -t 1-$NFILES
		cd \$SGE_O_WORKDIR

		module load krthornt/thorntonlab
		module load krthornt/SweepFinder/dev

		#infile=\`find /share/kevin2/krthornt/qtrait_paper/genome_scan -maxdepth 1 -name "*$sel.$timing"_shift"*.gz" |sort| head -n \$SGE_TASK_ID | tail -n 1\`
		infile=\`head -n \$SGE_TASK_ID /share/kevin2/krthornt/qtrait_paper/genome_scan/neutral_post.txt | tail -n 1\`
		stub=\`basename \$infile .gz\`
		outfile=\$stub.sweepfinder.txt
		if [ ! -e \$outfile.gz ]
		then
			mkdir \$stub
			cd \$stub
			NREPS=\`zgrep -c segs \$infile\`
			SweepFinder -m 100 <({ echo \"ms 100 \$NREPS -t 100\" & /share/kevin2/krthornt/qtrait_paper/src/trimFixations /share/kevin2/krthornt/qtrait_paper/genome_scan/\$infile ; }) \$outfile > /dev/null
			gzip \$outfile
			mv \$outfile.gz ..
			cd ..
			rm -rf \$stub
			NRECORDS=\`zcat \$outfile.gz|wc -l\`
			if [ \$NRECORDS -ne 20011 ]
			then
				#There was some sort of data loss,
				#so we kill the output file so
				#that we can resubmit the array
				rm -f \$outfile.gz
			fi
		fi
		EOF
		qsub -N sf.$timing.$sel $SCRIPTNAME
	fi
	done
done
