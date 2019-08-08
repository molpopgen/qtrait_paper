import glob
import gzip
import fwdpy11

for i in glob.glob('*.gz'):
    print(i)
    with gzip.open(i, 'rb') as f:
        pop = fwdpy11.DiploidPopulation.load_from_pickle_file(f)
        for j, k in zip(pop.fixation_times, pop.fixations):
            print(k.g, j, k.s)
