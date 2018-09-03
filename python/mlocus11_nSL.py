"""
This is a do-over of the nSL
calculations.  Here, we normalize based
on variants in the first and last window
of each locus and each replicate.

Note: seeding? Check how it was done.
"""
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
from collections import namedtuple
from libsequence.polytable import SimData
from libsequence.windows import Windows
from libsequence.summstats import nSLiHS
from libsequence.parallel import Scheduler
from fwdpy11.sampling import sample_separate

DataRecord = namedtuple(
    "DataRecord", ['generation', 'repid', 'locus', 'window', 'mean_nSLz'])

TempRecord = namedtuple("TempRecord", ['locus', 'window', 'values'])


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


def get_summstats(pop, repid, nsam):
    """
    The 'genome scan' bit.
    """
    ind = np.random.choice(pop.N, nsam, replace=False)
    s = sample_separate(pop, ind)
    locus = 0

    # The procedure here is:
    # 1. Get raw nSL for all windows at all loci
    # 2. For each locus, collect values from first/last window
    #    to use for normalizing.
    # 3. For all other windows per locus, get mean z-score

    reference_daf = []
    reference_values = []
    temp = []
    for si, bi, locus_index in zip(s, pop.locus_boundaries, range(len(s))):
        neut, sel = si
        neut.extend([i for i in sel])
        neut = sorted(neut, key=lambda x: x[0])
        sd = SimData(neut)
        assert sd.size() == 2 * args.nsam, "sample size error"
        w = Windows(sd, 1.0, 1.0, bi[0], bi[1])
        for i in range(len(w)):
            raw_nSL = nSLiHS(w[i])
            # Filter out non-finite values
            # and values where derived allele
            # present fewer than 3 times.
            raw_nSL = [i for i in raw_nSL if np.isfinite(
                i[0]) == 1 and i[2] > 3]
            nSL = np.array(raw_nSL)
            if len(nSL) > 0:
                if i == 0 or i == len(w) - 1:
                    reference_values.extend(nSL[:, 0].tolist())
                    reference_daf.extend(nSL[:, 2].tolist())
                else:
                    temp.append(TempRecord(locus, i, nSL))
        locus += 1

    # bin the reference data
    rdaf = np.array(reference_daf)
    rdaf_bins = np.digitize(rdaf, np.arange(0, 2 * args.nsam, 10))
    rstats = np.array(reference_values)

    mean_sd = {}
    for b in set(rdaf_bins):
        w = np.where(rdaf_bins == b)[0]
        if len(w) > 0:
            m = rstats[w].mean()
            sdev = rstats[w].std()
            if np.isfinite(sdev) == 1 and sdev > 0.0:
                mean_sd[b] = (m, sdev)

    rv = []
    # package up the data
    for t in temp:
        tb = np.digitize(t.values[:, 2],
                         np.arange(0, 2 * args.nsam, 10))
        zscores_win = np.array([])
        for b in set(tb):
            w = np.where(tb == b)[0]
            if b in mean_sd:
                m = mean_sd[b][0]
                sdev = mean_sd[b][1]
                zscores = (t.values[:, 0][w] - m) / sdev
                zscores_win = np.concatenate((zscores_win, zscores))
        mz = zscores_win.mean()
        rv.append(DataRecord(pop.generation, repid, t.locus,
                             t.window, mz))
    return rv


def process_replicate(argtuple):
    seed, args, infile = argtuple
    np.random.seed(seed)
    sched = Scheduler(1)
    tf = tarfile.open(args.tarfile, 'r')
    ti = tf.getmember(infile)
    lzma_file = tf.extract(ti)
    sstats = []
    with lzma.open(infile, 'rb') as f:
        while True:
            try:
                rec = pickle.load(f)
                g = rec[1].generation
                t1 = 9 * rec[1].N
                t2 = 14 * rec[1].N
                if g >= t1 and g <= t2:
                    ss = get_summstats(rec[1], rec[0], args.nsam)
                    sstats.extend(ss)
            except EOFError:
                print("returning {}".format(len(sstats)))
                os.remove(infile)
                return sstats
            except:
                print("unexpected exception")
                os.remove(infile)
                raise


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
    raw_args = [(i, args, j) for i, j in zip(
        sorted(np.random.choice(int(4e6), len(files), replace=False)), files)]
    # x = process_replicate(raw_args[0])
    # df = pd.DataFrame(x, columns=DataRecord._fields)
    # print(df.head())
    # sys.exit(0)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_replicate, i): i for i in raw_args}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            df = pd.DataFrame(fn, columns=DataRecord._fields)
            with sqlite3.connect(args.outfile) as conn:
                df.to_sql('data', conn, if_exists='append', index=False)
