# Calculate the fraction of large-
# and small- effect mutations, with
# respect to trait value.
import sys
import os
import pickle
import tarfile
import argparse
import lzma
import gzip
import math
import concurrent.futures
import numpy as np
import pandas as pd
import sqlite3


def make_parser():
    parser = argparse.ArgumentParser(
        description="Process multi-locus outputs")

    parser.add_argument('--tarfile', '-t', type=str,
                        help=".tar file name containing replicates")
    parser.add_argument('--nprocs', '-p', type=int, default=1,
                        help="Number of processes to use")
    parser.add_argument('--outfile', '-o', type=str, default=None,help="Output file name.")
    parser.add_argument('--mu', '-m', type=float, help="Mutation rate")
    parser.add_argument('--opt', '-O', type=float, help="New optimum value")
    return parser


def process_replicate(argtuple):
    """
    """
    args, infile = argtuple
    tf = tarfile.open(args.tarfile, 'r')
    ti = tf.getmember(infile)
    lzma_file = tf.extract(ti)
    optimum = 0.0
    rv = []

    # Mutations with effect sizes > threshold
    # are large-effect. Else, they are small.
    # Nuance is needed to note that this should
    # refer to |esize| of mutation.
    threshold = 2.0 * math.sqrt(2.0) * math.sqrt(args.mu)
    with lzma.open(infile, 'rb') as f:
        while True:
            try:
                repid, pop = pickle.load(f)
                if pop.generation < 8 * pop.N:
                    continue
                if pop.generation >= 10 * pop.N:
                    optimum = args.opt
                traits = np.array(pop.diploids.trait_array())
                mean_pheno = traits['g'].mean()
                muts = np.array(pop.mutations.array())
                muts_df = pd.DataFrame(muts)
                muts_df.loc[:, 'mcounts'] = pop.mcounts
                # Reduce to extant variants affecting trait value:
                muts_df = muts_df.query('mcounts > 0 and neutral == 0')
                # Mark variants as large-effect or not:
                muts_df.loc[:, 'large'] = np.abs(muts_df.s) > threshold
                # Mark mutations that just arose:
                muts_df.loc[:, 'new_mutation'] = (muts_df.g == pop.generation)

                # Slice up the DF a bunch of ways:
                new_large = muts_df.query('large == 1 and new_mutation == 1')
                new_small = muts_df.query('large == 0 and new_mutation == 1')
                old_large = muts_df.query('large == 1 and new_mutation == 0')
                old_small = muts_df.query('large == 0 and new_mutation == 0')

                # Distance pop still has to adapt
                # If distance is < 0, the population
                # has "over-shot" the optimum.
                distance = optimum - mean_pheno

                crit_value = abs(distance) / 2.0
                if distance > 0:
                    q = 's > 0 and s <= {}'.format(crit_value)
                else:
                    q = 's < 0 and s <= {}'.format(crit_value)
                if threshold > crit_value:
                    new_large_sweetspot = 0
                else:
                    new_large_sweetspot = len(new_large.query(q))

                new_small_sweetspot = len(new_small.query(q))
                # Add a summary as a big dict:
                rv.append({'rep': repid,
                           'generation': pop.generation,
                           'new_large': len(new_large),
                           'new_small': len(new_small),
                           'old_large': len(old_large),
                           'old_small': len(old_small),
                           'new_large_mean_esize': new_large.s.mean(),
                           'new_small_mean_esize': new_small.s.mean(),
                           'old_large_mean_esize': old_large.s.mean(),
                           'old_small_mean_esize': old_small.s.mean(),
                           'old_large_mean_freq': float(old_large.mcounts.mean()) / float(2 * pop.N),
                           'old_small_mean_freq': float(old_small.mcounts.mean()) / float(2 * pop.N),
                           'new_large_sweetspot': new_large_sweetspot,
                           'new_small_sweetspot': new_small_sweetspot
                           })
            except:
                break
    os.remove(infile)
    return (rv)


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    if os.path.exists(args.tarfile) is False:
        raise RuntimeError(args.tarfile + " could not be found")

    tf = tarfile.open(args.tarfile, 'r')
    files = tf.getnames()
    tf.close()
    raw_args = [(args, j) for j in files]
    of = gzip.open(args.outfile,'wb')
    of.close()
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.nprocs) as executor:
        futures = {executor.submit(process_replicate, i)                   : i for i in raw_args}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            fndf = pd.DataFrame(fn)
            fndf['mu'] = args.mu
            fndf['opt'] = args.opt
            fndf.to_csv(args.outfile,compression='gzip',mode='a',sep='\t',index=False)
