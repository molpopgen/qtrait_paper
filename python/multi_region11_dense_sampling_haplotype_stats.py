# 1. 10 loci separated by 50cM
# 2. Each locus is 10*theta and 10*rho long, where theta & rho are per-locus input params
# 3. The middle of each locus is where causative mutations happen at rate mu with
#    and sigmamu=0.25 and VS=1
# This script samples each generation for a while in order to capture *very(
# short-term dynamics of nSL, etc/

import fwdpy11 as fp11
import fwdpy11.wright_fisher
import fwdpy11.model_params
import fwdpy11.genetic_values
import fwdpy11.multilocus as fp11ml
import libsequence.polytable
import libsequence.summstats
from collections import namedtuple
import numpy as np
import pandas as pd
import sys
import datetime
import warnings
import argparse
import random
import pickle
import os
import lzma
import tarfile
import concurrent.futures
#import multiprocessing as mp

assert fwdpy11.__version__ > '0.1.4', "fwdpy11 version too old."


def make_parser():
    parser = argparse.ArgumentParser(
        description="Simple multi-locus simulation of single optimum shift")

    parser.add_argument('--N', '-N', type=int, default=5000,
                        help="Diploid population size")
    parser.add_argument('--mu', '-m', type=float, default=1e-3,
                        help="mutation rate to mutations affecting trait value")
    parser.add_argument('--theta', '-t', type=float,
                        default=100.0, metavar="THETA", help="4Nu")
    parser.add_argument('--rho', '-r', type=float,
                        default=100.0, metavar="RHO", help="4Nr")
    parser.add_argument('--seed', '-S', type=int,
                        default=None, help="RNG seed")
    parser.add_argument('--ncores', '-n', type=int, default=1)
    parser.add_argument('--nreps', '-nr', type=int, default=1)
    parser.add_argument('--nloci', '-nl', type=int, default=10)
    parser.add_argument('--sigmu', '-si', type=float, default=0.25)
    parser.add_argument('--opt', '-o', type=float, default=0.0)
    parser.add_argument('--outfile', type=str, default=None)
    return parser


class HapStatSampler(object):
    def __init__(self, repid, N):
        self.repid = repid
        self.N = N
        self.data = []

    def __call__(self, pop):
        c1 = (pop.generation >= 8 * self.N and pop.generation <
              10 * self.N and pop.generation % 500 == 0.0)
        c1b = (pop.generation >= 10 * self.N and
               pop.generation <= 10 * self.N + 500)
        c2 = (pop.generation >= 10 * self.N and pop.generation <=
              12 * self.N and pop.generation % 50 == 0.0)
        c3 = (pop.generation > 12 * self.N and pop.generation <=
              14 * self.N and pop.generation % 250 == 0.0)
        c4 = (pop.generation > 14 * self.N and pop.generation % 500 == 0.0)
        if c1 or c1b or c2 or c3 or c4:
            pass


def run_replicate(argtuple):
    args, repid, repseed = argtuple
    rng = fp11.GSLrng(repseed)
    NANC = args.N
    locus_boundaries = [(float(i + i * 11), float(i + i * 11 + 11))
                        for i in range(args.nloci)]
    nregions = [[fp11.Region(j[0], j[1], args.theta / (4. * float(NANC)), coupled=True)]
                for i, j in zip(range(args.nloci), locus_boundaries)]
    recregions = [[fp11.Region(j[0], j[1], args.rho / (4. * float(NANC)), coupled=True)]
                  for i, j in zip(range(args.nloci), locus_boundaries)]
    sregions = [[fp11.GaussianS(j[0] + 5., j[0] + 6., args.mu, args.sigmu, coupled=False)]
                for i, j in zip(range(args.nloci), locus_boundaries)]
    interlocus_rec = fp11ml.binomial_rec(rng, [0.5] * (args.nloci - 1))
    nlist = np.array([NANC] * 20 * NANC, dtype=np.uint32)
    env = [(0, 0, 1), (10 * NANC, args.opt, 1)]
    gv = fwdpy11.genetic_values.MlocusAdditive(
        2.0, fwdpy11.genetic_values, GSSmo(env))
    mutrates_s = [args.mu / float(args.nloci)] * args.nloci
    mutrates_n = [10 * args.theta / float(4 * NANC)] * args.nloci
    recrates = [10 * args.rho / float(4 * NANC)] * args.nloci
    pdict = {'nregions': nregions,
             'sregions': sregions,
             'recregions': recregions,
             'demography': nlist,
             'interlocus': interlocus_rec,
             'rates': (mutrates_n, mutrates_s, recrates),
             'gvalue': gv,
             'prune_selected': False
             }
    params = fp11.model_params.MlocusParamsQ(**pdict)
    recorder = HapStatSampler(repid, NANC)
    pop = fp11.MlocusPop(NANC, locus_boundaries)
    assert pop.nloci == len(locus_boundaries), "nloci != len(locus_boundaries)"
    fwdpy11.wright_fisher.evolve(rng, pop, params, recorder)
    # Make sure last gen got pickled!
    return repid, recorder.data


if __name__ == "__main__":
    # TODO list
    # 1. Track fixations to a file
    # 2. Track nSL and Garud stats to different file
    # 3. Or, do the analysis ahead of time! :)
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    initial_seed = args.seed
    np.random.seed(initial_seed)
    repseeds = np.random.randint(0, 42000000, args.nreps)
    arglist = [(args, i, j) for i, j in zip(range(len(repseeds)), repseeds)]
    tfn = args.tarfile
    if tfn is not None:
        tf = tarfile.open(tfn, "w")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.ncores) as executor:
        futures = {executor.submit(run_replicate, i): i for i in arglist}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            if tfn is not None:
                tf.add(fn)
                os.remove(fn)
    if tfn is not None:
        tf.close()
    # run_replicate(arglist[0])
    #P = mp.Pool(args.ncores)
    #res = P.imap_unordered(run_replicate,arglist)
    # P.close()
    # P.join()
