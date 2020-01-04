import gzip
import sys
import fwdpy11

with gzip.open(sys.argv[1],'rb') as f:
    pop = fwdpy11.DiploidPopulation.load_from_pickle_file(f)
