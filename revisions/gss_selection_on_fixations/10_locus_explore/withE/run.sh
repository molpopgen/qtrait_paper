#!/bin/bash

parallel python3 simulate_10loci_withE.py --popsize 5000 --VS 1 --mu 0.005 --filename himu.{}.gz --seed {} :::: himu &
parallel python3 simulate_10loci_withE.py --popsize 5000 --VS 1 --mu 0.001 --filename midmu.{}.gz --seed {} :::: midmu &
parallel python3 simulate_10loci_withE.py --popsize 5000 --VS 1 --mu 0.00025 --filename lowmu.{}.gz --seed {} :::: lowmu

ls -1 himu.*.gz midmu.*.gz lowmu.*.gz > infiles
cat infiles | sed 's/gz/sqlite3/' > outfiles
parallel --xapply python3 ../process_population.py --infile {1} --outfile {2} :::: infiles :::: outfiles
ls -1 himu.*.gz > himufiles.txt
ls -1 midmu.*.gz > midmufiles.txt
ls -1 lowmu.*.gz > lowmufiles.txt

for i in himufiles.txt midmufiles.txt lowmufiles.txt
do
    n=`basename $i .txt`
    cat $i | sed 's/gz/sqlite3/' > $n.outfiles.txt
    parallel --xapply python3 ../process_population.py --infile {1} --outfile {2} :::: $i :::: $n.outfiles.txt
done
