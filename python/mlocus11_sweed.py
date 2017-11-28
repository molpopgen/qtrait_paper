import sys
import os
import pickle
import tarfile
import argparse
import lzma
import concurrent.futures
import numpy as np
import pandas as pd
import sqlite3
import os
from libsequence.polytable import SimData
from fwdpy11.sampling import sample_separate
import collections
import subprocess
import shutil

def make_parser():
    parser = argparse.ArgumentParser(
        description="Run Sweed on each locus at each pickled time point.")

    parser.add_argument('--tarfile', '-t', type=str,
                        help=".tar file name containing replicates")
    parser.add_argument('--nprocs', '-p', type=int, default=1,
                        help="Number of processes to use")
    parser.add_argument('--outfile', '-o', type=str, help="Output file name.")
    parser.add_argument('--nsam', '-n', type=int, default=50,
                        help="Sample size (no. diploids)")
    parser.add_argument('--seed', '-s', type=int,
                        default=42, help="Random number seed")
    parser.add_argument('--opt',type=float,help="Optimum")
    parser.add_argument('--mu',type=float,help="Mutation rate")
    parser.add_argument('--gridsize',type=int,help="Grid size for SweeD",default=100)
    return parser


DataPoint = collections.namedtuple(
    'DataPoint', ['repid','generation', 'locus', 'LR', 'position', 'alpha'])
SweedEntry = collections.namedtuple(
    'SweedEntry', ['position', 'x', 'n', 'folded'])


def process_replicate(argtuple):
    seed, tdir, args, infile = argtuple
    np.random.seed(seed)
    if os.path.exists(tdir) is False:
        os.mkdir(tdir)
    dataFile = 'in.txt'
    #resultsFile =  'out.txt'
    resultsFile = 'SweeD_Report.run'
    tf = tarfile.open(args.tarfile, 'r')
    ti = tf.getmember(infile)
    lzma_file = tf.extract(ti)
    rv = []
    with lzma.open(ti.name, 'rb') as picklefile:
        while True:
            try:
                repid, pop = pickle.load(picklefile)
                ind = np.random.choice(pop.N, args.nsam, replace=False)
                s = pop.sample_ind(ind,separate=False)
                for locus, si in zip(range(len(s)), s):
                    sweed_data = []
                    for di in si:
                        sweed_data.append(SweedEntry(
                            di[0] * 1e6, di[1].count("1"), len(di[1]), 0))
                    data_df = pd.DataFrame(sweed_data)
                    with open(tdir + '/' + dataFile, 'w') as f:
                        data_df.to_csv(f, sep='\t', index=False)

                    x = subprocess.run(['SweeD', '-name', 'run', '-input', dataFile,
                                    '-grid', str(args.gridsize)], stdout=subprocess.PIPE, cwd=tdir)

                    res = pd.read_csv(tdir + '/' + resultsFile,
                                      sep='\t', skiprows=2)

                    LRmax = res['Likelihood'].argmax()
                    LR = res['Likelihood'].iloc[LRmax]
                    alpha = res['Alpha'].iloc[LRmax]
                    pos = res['Position'].iloc[LRmax]/1e6
                    rv.append(DataPoint(repid,pop.generation,locus,LR,pos,alpha))
                    os.remove(tdir + '/' + dataFile)
                    os.remove(tdir + '/' + resultsFile)
            except:
                break
    os.remove(infile)
    shutil.rmtree(tdir)
    return rv


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    np.random.seed(args.seed)
    if os.path.exists(args.tarfile) is False:
        raise RuntimeError(args.tarfile + " could not be found")

    if os.path.exists(args.outfile):
        os.remove(args.outfile)

    tf = tarfile.open(args.tarfile, 'r')
    files = tf.getnames()
    tf.close()
    raw_args = [(i, 'tdir' + str(args.mu) + '.' + str(args.opt) + '.' + str(i), args, j) for i, j in zip(
        sorted(np.random.choice(int(4e6), len(files), replace=False)), files)]
    conn = sqlite3.connect(args.outfile)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.nprocs) as executor:
        futures = {executor.submit(process_replicate, i): i for i in raw_args}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            data = pd.DataFrame(fn)
            data['opt'] = [args.opt]*len(data.index)
            data['mu'] = [args.mu]*len(data.index)

            data.to_sql('data',conn,index=False,if_exists='append')
    conn.close()
