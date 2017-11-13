import os
import shutil
import collections
import msprime
import pandas as pd
import subprocess
import numpy as np
import concurrent.futures
import sqlite3
import sys


SEED = 101
THETA = 1000.
RHO = 1000.
NSAM = 100
OUTFILE = ('sweed_nulldist.db')

SweedResult = collections.namedtuple('SweedResult', ['LR', 'position', 'alpha'])
SweedEntry = collections.namedtuple(
    'SweedEntry', ['position', 'x', 'n', 'folded'])


def do_work(args):
    seed, nreps, tdir = args
    if os.path.exists(tdir) is False:
        os.mkdir(tdir)

    dataFile = 'in.txt'
    resultsFile = 'SweeD_Report.run'
    nulldist = []
    for ts in msprime.simulate(NSAM, mutation_rate=THETA / 4.,
                               recombination_rate=RHO / 4.,
                               random_seed=int(seed), num_replicates=nreps):
        data = []
        for variant in ts.variants():
            data.append(SweedEntry(
                variant.position * 1e6, len(np.where(variant.genotypes == 1)[0]), NSAM, 0))
        data_df = pd.DataFrame(data)
        with open(tdir + '/' + dataFile, 'w') as f:
            data_df.to_csv(f, sep='\t', index=False)

        x = subprocess.run(['SweeD', '-name', 'run', '-input', dataFile,
                            '-grid', '100'], stdout=subprocess.PIPE, cwd=tdir)

        res = pd.read_csv(tdir + '/' + resultsFile,
                          sep='\t', skiprows=2)
        LRmax = res['Likelihood'].argmax()
        LR = res['Likelihood'].iloc[LRmax]
        alpha = res['Alpha'].iloc[LRmax]
        pos = res['Position'].iloc[LRmax] / 1e6
        nulldist.append(SweedResult(LR, pos, alpha))
    shutil.rmtree(tdir)
    return nulldist


if __name__ == "__main__":
    np.random.seed(SEED)
    seeds = [int(i) for i in np.random.choice(40000000, 72, replace=False)]
    nreps = [int(1e6 / 72.)] * 72
    nreps[0] = 1000000 - sum(nreps)
    if os.path.exists(OUTFILE):
        os.remove(OUTFILE)

    conn = sqlite3.connect(OUTFILE)
    with concurrent.futures.ProcessPoolExecutor(72) as executor:
        futures = {executor.submit(
            do_work, (i, j, "sweed_temp" + str(i))): i for i, j in zip(seeds, nreps)}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            data = pd.DataFrame(fn)

            data.to_sql('data', conn, index=False, if_exists='append')
    conn.close()
