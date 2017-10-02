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
from libsequence.summstats import PolySIM

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


def process_replicate(argtuple):
    seed, args, infile = argtuple
    np.random.seed(seed)
    tf = tarfile.open(args.tarfile, 'r')
    ti = tf.getmember(infile)
    lzma_file = tf.extract(ti)
    rv = []
    statnames = ['thetapi', 'tajimasd', 'hprime']
    with lzma.open(infile, 'rb') as f:
        while True:
            try:
                rep, pop = pickle.load(f)
                ind = np.random.choice(pop.N, args.nsam, replace=False)
                s = sample_separate(pop, ind)
                # Only loci 0, 5, and len(s)-1 can have
                # any neutral variants
                for locus in [0, 5, len(s) - 1]:
                    sd = SimData(s[locus][0])
                    ps = PolySIM(sd)
                    for stat in statnames:
                        d = {'rep': rep,
                             'locus': locus,
                             'generation': pop.generation,
                             'stat': stat,
                             'value': getattr(ps, stat)()
                             }
                        rv.append(d)
            except:
                break
    os.remove(infile)
    return rv


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    np.random.seed(args.seed)
    if os.path.exists(args.tarfile) is False:
        raise RuntimeError(args.tarfile + " could not be found")

    dbname = args.tarfile.replace('.tar', '.genome_scan.db')
    #         args.tarfile.replace('.tar','.genetic_values_per_locus.db'),
    #         args.tarfile.replace('.tar','.vg_over_time.db')]
    # for d in dbnames:
    #    if os.path.exists(d):
    #        os.remove(d)
    tf = tarfile.open(args.tarfile, 'r')
    files = tf.getnames()
    tf.close()
    raw_args = [(i, args, j) for i, j in zip(
        sorted(np.random.choice(int(4e6), len(files), replace=False)), files)]
    # x=process_replicate(raw_args[0])
    # print(x)
    # sys.exit(0)
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        futures = {executor.submit(process_replicate, i): i for i in raw_args}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            conn = sqlite3.connect(dbname)
            df = pd.DataFrame(fn)
            df.to_sql('stats', conn, if_exists="append", index=False)
            conn.close()
