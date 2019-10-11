#!/bin/bash

parallel --xapply python3 ../process_10_loci.py {1} {2} {3} :::: infiles_processing.txt :::: outfiles_processing.txt :::: seed_processing.txt
