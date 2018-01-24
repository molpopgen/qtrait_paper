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
import collections
from libsequence.polytable import SimData
from libsequence.windows import Windows
from libsequence.summstats import  nSLiHS
from libsequence.parallel import Scheduler
from fwdpy11.sampling import sample_separate


def make_parser():
    parser = argparse.ArgumentParser(
        description="Process multi-locus outputs")

    parser.add_argument('--tarfile', '-t', type=str,
                        help=".tar file name containing replicates")
    parser.add_argument('--nprocs', '-p', type=int, default=1,
                        help="Number of processes to use")
    parser.add_argument('--outfile', '-o', type=str, help="Output file name.")
    parser.add_argument('--nsam', '-n', type=int, default=50,
                        help="Sample size (no. diploids)")
    parser.add_argument('--seed', '-s', type=int,
                        default=42, help="Random number seed")
    return parser

Datum=collections.namedtuple("Datum",['repid','generation','locus','window','score','daf'])

def get_raw_nSL(pop, repid, nsam):
    """
    The 'genome scan' bit.
    """
    ind = np.random.choice(pop.N, nsam, replace=False)
    s = sample_separate(pop, ind)
    rv=[]
    for si, bi, locus_index in zip(s, pop.locus_boundaries, range(len(s))):
        neut, sel = si
        neut.extend([i for i in sel])
        neut = sorted(neut,key = lambda x: x[0])
        sd = SimData(neut)
        w = Windows(sd, 1.0, 1.0, bi[0], bi[1])
        for i in range(len(w)):
            ps = PolySIM(w[i])
            gs = garudStats(w[i])

            raw_nSL = nSLiHS(w[i])
            raw_nSL = [i for i in raw_nSL if np.isfinite(i[0]) == 1]
            x =[Datum(repid,pop.generation,locus_index,window,i[0],i[2]) for i in raw_nSL]
            rv.extend(x)
    return rv


def process_replicate(argtuple):
    seed, args, infile = argtuple
    np.random.seed(seed)
    sched = Scheduler(1)
    tf = tarfile.open(args.tarfile, 'r')
    ti = tf.getmember(infile)
    lzma_file = tf.extract(ti)
    nSL=[]
    with lzma.open(infile, 'rb') as f:
        while True:
            try:
                rec = pickle.load(f)
                nSL.extend(get_raw_nSL(pop,repid,args.nsam))
            except EOFError:
                os.remove(infile)
                return (nSL,)
            except:
                os.remove(infile)
                raise


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    np.random.seed(args.seed)
    if os.path.exists(args.tarfile) is False:
        raise RuntimeError(args.tarfile + " could not be found")

    dbnames = [args.tarfile.replace('.tar', '.nSL.db')]
    for d in dbnames:
        if os.path.exists(d):
            os.remove(d)
    tf = tarfile.open(args.tarfile, 'r')
    files = tf.getnames()
    tf.close()
    raw_args = [(i, args, j) for i, j in zip(
        sorted(np.random.choice(int(4e6), len(files), replace=False)), files)]
    # x=process_replicate(raw_args[0])
    # print(x)
    # sys.exit(0)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_replicate, i): i for i in raw_args}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            for dbname, data in zip(dbnames, fn):
                conn = sqlite3.connect(dbname)
                df = pd.DataFrame(data)
                df.to_sql('data', conn, if_exists="append", index=False)
                conn.close()
    for dbname in dbnames:
        conn = sqlite3.connect(dbname)
        conn.execute("create index gen_rep if not exists on data")
        conn.close()

