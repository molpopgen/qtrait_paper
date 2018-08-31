#!/bin/bash

for i in pops.*.tar
do
    n=`basename $i .tar`
    PYTHONPATH=$HOME/src/fwdpy11 python3 ../python/mlocus11_nSL.py -t $i  --nprocs 40 -o $n.nSLz.sqlite3
done
