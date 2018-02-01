# Simulate the null distribution of stats for a window
import msprime
import libsequence.msprime
import libsequence.summstats
import libsequence.parallel
import pandas as pd
import numpy as np
import sqlite3
from collections import namedtuple
import concurrent.futures
import os


# n is 100 b/c we sample 50 diploids
NSAM = 100
SEED = 42
THETA = 1000. / 11.
RHO = 1000. / 11.
NREPS = 50000  # per process

StatHolder = namedtuple(
    'StatHolder', ['pbig2', 'pbig3'])


def do_work(local_seed):
    nulldist = []
    sched = libsequence.parallel.Scheduler(1)
    for ts in msprime.simulate(NSAM, mutation_rate=THETA / 4.,
                               recombination_rate=RHO / 4.,
                               random_seed=int(local_seed), num_replicates=NREPS):
        d = libsequence.msprime.make_SimData(ts)
        raw_nSL = libsequence.summstats.nSLiHS(d)
        raw_nSL = [i for i in raw_nSL if np.isfinite(
            i[0]) == 1 and min(i[2], NSAM - i[2]) > 3]
        nSL = np.array(raw_nSL)
        bins = np.digitize(nSL[:, 2], np.arange(0, NSAM, 5))
        nbig2 = 0
        nbig3 = 0
        nscores=0
        for b in set(bins):
            binscores = nSL[:, 0][np.where(bins == b)[0]]
            if len(binscores) > 1:
                bmean = binscores.mean()
                sd = binscores.std()
                if sd > 0.0 and np.isfinite(sd):
                    binscores = (binscores - bmean) / sd
                    nscores+=len(binscores)
                    gt2 = np.where(np.abs(binscores >= 2.0))[0]
                    nbig2 += len(gt2)
                    gt3 = np.where(np.abs(binscores >= 3.0))[0]
                    nbig3 += len(gt3)
        nulldist.append(StatHolder(nbig2/nscores, nbig3/nscores))
    return nulldist


if __name__ == "__main__":
    np.random.seed(SEED)
    raw_args = [(i) for i in np.random.choice(1000000, 20, replace=False)]
    if os.path.exists("nulldist_abs_nSL.db"):
        os.remove("nulldist_abs_nSL.db")
    conn = sqlite3.connect("nulldist_abs_nSL.db")
    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(do_work, i): i for i in raw_args}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            data = pd.DataFrame(fn)
            data.to_sql('data', conn, index=False, if_exists='append')
    conn.close()
