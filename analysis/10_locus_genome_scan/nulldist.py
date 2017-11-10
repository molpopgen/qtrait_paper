# Simulate the null distribution of stats for a window
import msprime
import libsequence.msprime
import libsequence.summstats
import libsequence.parallel
import pandas as pd
import numpy as np
import sqlite3
from collections import namedtuple

# n is 100 b/c we sample 50 diploids
n = 100
seed = 42
theta = 1000. / 11.
rho = 1000. / 11.
nreps = 1000000

StatHolder = namedtuple(
    'StatHolder', ['tajd', 'hprime', 'H1', 'H12', 'H2H1', 'max_abs_nSL'])


nulldist = []
sched = libsequence.parallel.Scheduler(1)
for ts in msprime.simulate(n, mutation_rate=theta / 4.,
                           recombination_rate=rho / 4.,
                           random_seed=42, num_replicates=nreps):
    d = libsequence.msprime.make_SimData(ts)
    ad = libsequence.summstats.PolySIM(d)
    gs = libsequence.summstats.garudStats(d)
    raw_nSL = libsequence.summstats.nSLiHS(d)
    raw_nSL = [i for i in raw_nSL if np.isfinite(i[0]) == 1 and min(i[2],n-i[2]) > 3]
    nSL = np.array(raw_nSL)
    bins = np.digitize(nSL[:, 2], np.arange(0, n, 5))
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
    nulldist.append(StatHolder(ad.tajimasd(),ad.hprime(),gs['H1'],gs['H12'],gs['H2H1'],max_nSL))

df = pd.DataFrame(nulldist)
conn = seqlite3.connect("nulldist.db")
df.to_sql('data',conn,if_exists='replace')
