# Calculate the fraction of large-
# and small- effect mutations, with
# respect to trait value.
import sys
import os
import pickle
import tarfile
import argparse
import lzma
import math
import concurrent.futures
import numpy as np
import pandas as pd
import sqlite3
import os


def make_parser():
    parser = argparse.ArgumentParser(
        description="Process multi-locus outputs")

    parser.add_argument('--tarfile', '-t', type=str,
                        help=".tar file name containing replicates")
    parser.add_argument('--nprocs', '-p', type=int, default=1,
                        help="Number of processes to use")
    parser.add_argument('--outfile', '-o', type=str, help="Output file name.")
    parser.add_argument('--mu', '-m', type=float, help="Mutation rate")
    parser.add_argument('--opt', '-O', type=float, help="New optimum value")
    return parser


def process_replicate(argtuple):
    """
    """
    args, infile = argtuple
    np.random.seed(seed)
    sched = Scheduler(1)
    tf = tarfile.open(args.tarfile, 'r')
    ti = tf.getmember(infile)
    lzma_file = tf.extract(ti)
    optimum = 0.0
    rv = []

    # Mutations with effect sizes > threshold
    # are large-effect. Else, they are small.
    # Nuance is needed to note that this should
    # refer to |esize| of mutation.
    threshold = 2.0 * math.sqrt(2.0) * sqrt(args.mu)
    with lzma.open(infile, 'rb') as f:
        while True:
            try:
                repid, pop = pickle.load(f)
                if pop.generation >= 10 * pop.N:
                    optimum = args.opt
                traits = np.array(pop.diploids.trait_array())
                mean_pheno = traits['g'].mean()
                muts = np.array(pop.mutations.array())
                muts_df = pd.DataFrame(muts)
                muts_df['mcounts'] = pop.mcounts
                # Reduce to extant variants affecting trait value:
                muts_df = muts_df[(muts['mcounts'] > 0) &
                                  (muts.neutral is False)]
                # Mark variants as large-effect or not:
                muts_df['large'] = np.abs(muts_df.s) > threshold
                # Mark mutations that just arose:
                muts_df['new_mutation'] = (muts_df.g == pop.generation)

                # Distance pop still has to adapt
                # There is a problem due to over-shooting here!
                distance = optimum - mean_pheno

                # Slice up the DF a bunch of ways:
                new_large_positive = muts_df.query(
                    'large == 1 and new_mutation == 1 and s > 0')
                new_large_negative = muts_df.query(
                    'large == 1 and new_mutation == 1 and s < 0')
                new_small_positive = muts_df.query(
                    'large == 0 and new_mutation == 1 and s > 0')
                new_small_negative = muts_df.query(
                    'large == 0 and new_mutation == 1 and s < 0')

                old_large_positive = muts_df.query(
                    'large == 1 and new_mutation == 0 and s > 0')
                old_large_negative = muts_df.query(
                    'large == 1 and new_mutation == 0 and s < 0')
                old_small_positive = muts_df.query(
                    'large == 0 and new_mutation == 0 and s > 0')
                old_small_negative = muts_df.query(
                    'large == 0 and new_mutation == 0 and s < 0')

                # Add a summary as a big dict:
                rv.append({'rep': repid,
                           'generation': pop.generation,
                           'new_large_positive': len(new_large_positive),
                           'new_large_negative': len(new_large_negative),
                           'new_small_positive': len(new_small_positive),
                           'new_small_negative': len(new_small_negative),
                           'old_large_positive': len(old_large_positive),
                           'old_large_negative': len(old_large_negative),
                           'old_small_positive': len(old_small_positive),
                           'old_small_negative': len(old_small_negative),
                           'new_large_positive_mean_esize': new_large_positive.s.mean(),
                           'new_large_negative_mean_esize': new_large_negative.s.mean(),
                           'new_small_positive_mean_esize': new_small_positive.s.mean(),
                           'new_small_negative_mean_esize': new_small_negative.s.mean(),
                           'old_large_positive_mean_esize': old_large_positive.s.mean(),
                           'old_large_negative_mean_esize': old_large_negative.s.mean(),
                           'old_small_positive_mean_esize': old_small_positive.s.mean(),
                           'old_small_negative_mean_esize': old_small_negative.s.mean(),
                           'new_large_positive_mean_freq': new_large_positive.mcounts.mean() / float(2 * pop.N),
                           'new_large_negative_mean_freq': new_large_negative.mcounts.mean() / float(2 * pop.N),
                           'new_small_positive_mean_freq': new_small_positive.mcounts.mean() / float(2 * pop.N),
                           'new_small_negative_mean_freq': new_small_negative.mcounts.mean() / float(2 * pop.N),
                           'old_large_positive_mean_freq': old_large_positive.mcounts.mean() / float(2 * pop.N),
                           'old_large_negative_mean_freq': old_large_negative.mcounts.mean() / float(2 * pop.N),
                           'old_small_positive_mean_freq': old_small_positive.mcounts.mean() / float(2 * pop.N),
                           'old_small_negative_mean_freq': old_small_negative.mcounts.mean() / float(2 * pop.N)
                           })
                except:
                    break
    os.remove(infile)
    return (rv)


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    np.random.seed(args.seed)
    if os.path.exists(args.tarfile) is False:
        raise RuntimeError(args.tarfile + " could not be found")

    dbnames = []
    # for d in dbnames:
    #    if os.path.exists(d):
    #        os.remove(d)
    tf = tarfile.open(args.tarfile, 'r')
    files = tf.getnames()
    tf.close()
    raw_args = [(args, j) for j in files]
    # x=process_replicate(raw_args[0])
    # print(x)
    # sys.exit(0)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.nprocs) as executor:
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
