# 1. 10 loci separated by 50cM
# 2. Each locus is 10*theta and 10*rho long, where theta & rho are per-locus input params
# 3. The middle of each locus is where causative mutations happen at rate mu with
#   and sigmamu=0.25 and VS=1
# This script differs from multi_region11.py in that it doesn't pickle each pop.
# Rather, it records a few summary stats and moves on.
import fwdpy11 as fp11
import fwdpy11.wright_fisher_qtrait as fp11qt
import fwdpy11.model_params as fp11mp
import fwdpy11.trait_values as fp11tv
import fwdpy11.multilocus as fp11ml
from fwdpy11.sampling import sample_separate
from libsequence.polytable import SimData
from libsequence.windows import Windows
from libsequence.summstats import PolySIM, garudStats, nSLiHS
import numpy as np
import pandas as pd
import sys
import gc
import math
import datetime
import warnings
import argparse
import random
import pickle
import os
import sqlite3
import concurrent.futures
import scipy.stats
import scipy.optimize
from collections import namedtuple


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
    parser.add_argument('--plarge', type=float, default=None,
                        help="Probability of a large effect mutation. Calculated once at start, based on --mu and VS = 1")
    parser.add_argument('--vsopt', type=float, default=1.0,
                        help="VS after optimum shift")
    parser.add_argument('--nsam', type=int, default=50,
                        help="Sample size (no. diploids)")
    parser.add_argument('--statsdb', type=str, default=None,
                        help="File name for population stats database")
    parser.add_argument('--scandb', type=str, default=None,
                        help="File name for genome scan database")
    return parser


DataPoint = namedtuple('DataPoint', ['generation', 'VG', 'zbar', 'wbar'])
GenomeScanDataPoint = namedtuple(
    'GenomeScanDataPoint', ['generation', 'locus', 'window', 'D', 'Hprime'])


def get_summstats(pop, nsam, generation):
    """
    The 'genome scan' bit.
    """
    ind = np.random.choice(pop.N, nsam, replace=False)
    s = sample_separate(pop, ind)
    rv = []
    locus = 0
    for si, bi, locus_index in zip(s, pop.locus_boundaries, range(len(s))):
        neut, sel = si
        neut.extend([i for i in sel])
        neut = sorted(neut, key=lambda x: x[0])
        sd = SimData(neut)
        w = Windows(sd, 1.0, 1.0, bi[0], bi[1])
        for i in range(len(w)):
            ps = PolySIM(w[i])
            rv.append(GenomeScanDataPoint(
                generation, locus, i, ps.tajimasd(), ps.hprime()))
        locus += 1
    return rv


class Recorder(object):
    def __init__(self, repid, N, nsam):
        self.repid = repid
        self.N = N
        self.nsam = nsam
        self.data = []
        self.genome_scan_data = []

    def __call__(self, pop):
        c1 = (pop.generation >= 8 * self.N and pop.generation <
              10 * self.N and pop.generation % 500 == 0.0)
        c2 = (pop.generation >= 10 * self.N and pop.generation <=
              12 * self.N and pop.generation % 50 == 0.0)
        c3 = (pop.generation > 12 * self.N and pop.generation <=
              14 * self.N and pop.generation % 250 == 0.0)
        c4 = (pop.generation > 14 * self.N and pop.generation % 500 == 0.0)
        if c1 or c2 or c3 or c4:
            traits = np.array(pop.diploids.trait_array())
            self.data.append(DataPoint(
                pop.generation, traits['g'].var(), traits['g'].mean(), traits['w'].mean()))
            self.genome_scan_data.extend(
                get_summstats(pop, self.nsam, pop.generation))


def gamma_hat(VS, mu):
    return 2.0 * math.sqrt(2.0) * math.sqrt(VS * mu)


def generate_gaussian_function_to_minimize(ghat, plarge):
    return lambda x: abs(2.0 * (1.0 - scipy.stats.norm.cdf(ghat, scale=x)) - plarge)


def get_gaussian_sigma(F):
    res = scipy.optimize.minimize_scalar(F, bounds=(0, 100), method='bounded')
    return res.x


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

    # Get the variance in effect sizes
    sigmu = args.sigmu  # default
    if args.plarge is not None:
        ghat = gamma_hat(1.0, args.mu)
        F = generate_gaussian_function_to_minimize(ghat, args.plarge)
        sigmu = get_gaussian_sigma(F)

    sregions = [[fp11.GaussianS(j[0] + 5., j[0] + 6., args.mu, sigmu, coupled=False)]
                for i, j in zip(range(args.nloci), locus_boundaries)]
    interlocus_rec = fp11ml.binomial_rec(rng, [0.5] * (args.nloci - 1))
    nlist = np.array([NANC] * 20 * NANC, dtype=np.uint32)
    env = [(0, 0, 1), (10 * NANC, args.opt, args.vsopt)]
    pdict = {'nregions': nregions,
             'sregions': sregions,
             'recregions': recregions,
             'demography': nlist,
             'interlocus': interlocus_rec,
             'agg': fp11ml.AggAddTrait(),
             'gvalue': fp11ml.MultiLocusGeneticValue([fp11tv.SlocusAdditiveTrait(2.0)] * args.nloci),
             'trait2w': fp11qt.GSSmo(env),
             'mutrates_s': [args.mu / float(args.nloci)] * args.nloci,
             'mutrates_n': [10 * args.theta / float(4 * NANC)] * args.nloci,
             'recrates': [10 * args.rho / float(4 * NANC)] * args.nloci,
             'prune_selected': False
             }
    params = fp11.model_params.MlocusParamsQ(**pdict)
    recorder = Recorder(repid, NANC, args.nsam)
    pop = fp11.MlocusPop(NANC, args.nloci, locus_boundaries)
    fp11qt.evolve(rng, pop, params, recorder)
    return recorder


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    initial_seed = args.seed
    np.random.seed(initial_seed)
    repseeds = np.random.randint(0, 42000000, args.nreps)
    arglist = [(args, i, j) for i, j in zip(range(len(repseeds)), repseeds)]
    if os.path.exists(args.statsdb):
        os.remove(args.statsdb)
    if os.path.exists(args.scandb):
        os.remove(args.scandb)

    statscon = sqlite3.connect(args.statsdb)
    scancon = sqlite3.connect(args.scandb)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.ncores) as executor:
        futures = {executor.submit(run_replicate, i): i for i in arglist}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            vgstats = pd.DataFrame(fn.data)
            vgstats['repid'] = [fn.repid] * len(vgstats.index)
            vgstats.to_sql('data', statscon, if_exists='append', index=False)
            scanstats = pd.DataFrame(fn.genome_scan_data)
            scanstats['repid'] = [fn.repid] * len(scanstats.index)
            scanstats.to_sql('data', scancon, if_exists='append', index=False)
    statscon.close()
    scancon.close()
