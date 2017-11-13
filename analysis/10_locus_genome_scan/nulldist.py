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
    'StatHolder', ['tajd', 'hprime', 'H1', 'H12', 'H2H1', 'max_abs_nSL'])


def do_work(local_seed):
    nulldist = []
    sched = libsequence.parallel.Scheduler(1)
    for ts in msprime.simulate(NSAM, mutation_rate=THETA / 4.,
                               recombination_rate=RHO / 4.,
                               random_seed=int(local_seed), num_replicates=NREPS):
        d = libsequence.msprime.make_SimData(ts)
        ad = libsequence.summstats.PolySIM(d)
        gs = libsequence.summstats.garudStats(d)
        raw_nSL = libsequence.summstats.nSLiHS(d)
        raw_nSL = [i for i in raw_nSL if np.isfinite(
            i[0]) == 1 and min(i[2], NSAM - i[2]) > 3]
        nSL = np.array(raw_nSL)
        bins = np.digitize(nSL[:, 2], np.arange(0, NSAM, 5))
        max_nSL = np.nan
        for b in set(bins):
            binscores = nSL[:, 0][np.where(bins == b)[0]]
            if len(binscores) > 1:
                bmean = binscores.mean()
                sd = binscores.std()
                if sd > 0.0:
                    binscores = (binscores - bmean) / sd
                    bmax = binscores[np.where(
                        np.abs(binscores) == max(np.abs(binscores)))[0]][0]
                    if np.isnan(max_nSL) or np.abs(bmax) > np.abs(max_nSL):
                        max_nSL = bmax
        nulldist.append(StatHolder(ad.tajimasd(), ad.hprime(),
                                   gs['H1'], gs['H12'], gs['H2H1'], max_nSL))
    return nulldist


if __name__ == "__main__":
    np.random.seed(SEED)
    raw_args = [(i) for i in np.random.choice(1000000, 20, replace=False)]
    if os.path.exists("nulldist.db"):
        os.remove("nulldist.db")
    conn = sqlite3.connect("nulldist.db")
    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(do_work, i): i for i in raw_args}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            data = pd.DataFrame(fn)
            data.to_sql('data', conn, index=False, if_exists='append')
    conn.close()
