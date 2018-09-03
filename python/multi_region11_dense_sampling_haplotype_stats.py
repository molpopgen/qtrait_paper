# 1. 10 loci separated by 50cM
# 2. Each locus is 10*theta and 10*rho long, where theta & rho are per-locus input params
# 3. The middle of each locus is where causative mutations happen at rate mu with
#    and sigmamu=0.25 and VS=1
# This script samples each generation for a while in order to capture *very*
# short-term dynamics of nSL, etc.

import fwdpy11 as fp11
import fwdpy11.wright_fisher
import fwdpy11.model_params
import fwdpy11.genetic_values
import fwdpy11.sampling
import fwdpy11.multilocus as fp11ml
import libsequence.polytable
import libsequence.summstats
import libsequence.windows
from collections import namedtuple
import numpy as np
import pandas as pd
import sys
import argparse
import sqlite3
import os
import concurrent.futures
# import multiprocessing as mp

assert fwdpy11.__version__ > '0.1.4', "fwdpy11 version too old."

TempRecord = namedtuple("TempRecord", ['locus', 'window', 'values'])
DataRecord = namedtuple(
    "DataRecord", ['generation', 'repid', 'locus', 'window',
                   'mean_nSLz', 'H1', 'H12', 'H2H1'])


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
    parser.add_argument('--nsam', type=int, default=None,
                        help="Sample size (no. diploids)")
    parser.add_argument('--simlen', type=int, default=20,
                        help="Simulation length, in multiples of N generations (int).")
    parser.add_argument('--ncores', '-n', type=int, default=1)
    parser.add_argument('--nreps', '-nr', type=int, default=1)
    parser.add_argument('--nloci', '-nl', type=int, default=10)
    parser.add_argument('--sigmu', '-si', type=float, default=0.25)
    parser.add_argument('--opt', '-o', type=float, default=0.0)
    parser.add_argument('--outfile', type=str, default=None,
                        help="sqlite3 output file for summary stats")
    parser.add_argument('--fixationsfile', type=str,
                        default=None, help="sqlite3 file for fixation data")
    return parser


class HapStatSampler(object):
    def __init__(self, repid, N, nsam):
        self.repid = repid
        self.N = N
        self.nsam = nsam
        self.data = []

    def __call__(self, pop):
        c1 = (pop.generation >= 8 * self.N and pop.generation <
              10 * self.N and pop.generation % 500 == 0.0)
        # This is the new thing: we sample every generation
        # starting from the optimum shift to 0.1N gens later...
        c1b = (pop.generation >= 10 * self.N and
               pop.generation <= 10 * self.N + int(0.1*float(self.N)))
        # ...and then revert back to every 50 generations
        c2 = (pop.generation >= 10 * self.N and
              pop.generation <= 12 * self.N and pop.generation % 50 == 0.0)
        c3 = (pop.generation > 12 * self.N and pop.generation <=
              14 * self.N and pop.generation % 250 == 0.0)
        c4 = (pop.generation > 14 * self.N and pop.generation % 500 == 0.0)
        if c1 or c1b or c2 or c3 or c4:
            # TODO:
            # 1. Get the data matrix
            # 2. Convert it to a "sample"
            # 3. split by locus
            # 4. Create sliding windows
            # 5. get stats, copying the normaliation routing from mlocus11_nSL.py
            ind = np.random.choice(pop.N, self.nsam, replace=False)
            dm = pop.sample(ind, True)

            # Note: all of the below is a big performance hit
            # due to not having "libseq 2.0"!
            s = fwdpy11.sampling.matrix_to_sample(dm)
            neut_per_locus = fwdpy11.sampling.separate_samples_by_loci(
                pop.locus_boundaries, s[0])
            sel_per_locus = fwdpy11.sampling.separate_samples_by_loci(
                pop.locus_boundaries, s[1])
            # Store data for binning nSL from
            # first and last window, pooled
            # across each locus
            reference_daf = []
            reference_values = []
            # Store nSL values at all internal windows
            # for binning later
            temp = []
            # Intermediate structure for data records
            tempData = []
            for neut, sel, bi, locus in zip(neut_per_locus, sel_per_locus,
                                            pop.locus_boundaries, range(len(neut_per_locus))):
                neut.extend([i for i in sel])
                neut = sorted(neut, key=lambda x: x[0])
                sd = libsequence.polytable.SimData(neut)
                w = libsequence.windows.Windows(sd, 1.0, 1.0, bi[0], bi[1])

                # Iterate over each window
                for i in range(len(w)):
                    gs = libsequence.summstats.garudStats(w[i])
                    tempData.append(DataRecord(
                        pop.generation, self.repid, locus, i,
                        np.nan, gs['H1'], gs['H12'], gs['H2H1']))
                    raw_nSL = libsequence.summstats.nSLiHS(w[i])
                    # Filter out non-finite values
                    # and values where derived allele
                    # present fewer than 3 times.
                    raw_nSL = [i for i in raw_nSL if np.isfinite(
                        i[0]) == 1 and i[2] > 3]
                    if len(raw_nSL) > 0:
                        nSL = np.array(raw_nSL)
                        if i == 0 or i == len(w) - 1:
                            reference_values.extend(nSL[:, 0].tolist())
                            reference_daf.extend(nSL[:, 2].tolist())
                        else:
                            temp.append(TempRecord(locus, i, nSL))

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
                tempData[t.window] = tempData[t.window]._replace(mean_nSLz=mz)

            self.data.extend(tempData)


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
    interlocus_rec = fp11ml.binomial_rec([0.5] * (args.nloci - 1))
    nlist = np.array([NANC] * args.simlen * NANC, dtype=np.uint32)

    # the tuples are (time, z_o, VS)
    env = [(0, 0, 1), (10 * NANC, args.opt, 1)]
    gv = fwdpy11.genetic_values.MlocusAdditive(
        2.0, fwdpy11.genetic_values.GSSmo(env))
    mutrates_s = [args.mu / float(args.nloci)] * args.nloci
    mutrates_n = [10 * args.theta / float(4 * NANC)] * args.nloci
    recrates = [10 * args.rho / float(4 * NANC)] * args.nloci
    pdict = {'nregions': nregions,
             'sregions': sregions,
             'recregions': recregions,
             'demography': nlist,
             'interlocus_rec': interlocus_rec,
             'rates': (mutrates_n, mutrates_s, recrates),
             'gvalue': gv,
             'prune_selected': False
             }
    params = fp11.model_params.ModelParams(**pdict)
    recorder = HapStatSampler(repid, NANC, args.nsam)
    pop = fp11.MlocusPop(NANC, locus_boundaries)
    assert pop.nloci == len(locus_boundaries), "nloci != len(locus_boundaries)"
    fwdpy11.wright_fisher.evolve(rng, pop, params, recorder)
    return repid, recorder.data, pop


if __name__ == "__main__":
    # TODO list
    # 1. Track fixations to a file
    # 2. Track nSL and Garud stats to different file
    # 3. Or, do the analysis ahead of time! :)
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    assert args.nsam is not None, "nsam cannot be None"
    assert args.outfile is not None, "outfile cannot be None"
    assert args.fixationsfile is not None, "fixationsfile cannot be None"
    assert args.simlen > 10, "Simulation must be at least 10N generations long."

    for i in [args.outfile, args.fixationsfile]:
        if os.path.exists(i):
            os.remove(i)

    initial_seed = args.seed
    np.random.seed(initial_seed)
    repseeds = np.random.randint(0, 42000000, args.nreps)
    while len(np.unique(repseeds)) != args.nreps:
        repseeds = np.random.randint(0, 42000000, args.nreps)
    arglist = [(args, i, j) for i, j in zip(range(len(repseeds)), repseeds)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.ncores) as executor:
        futures = {executor.submit(run_replicate, i): i for i in arglist}
        for fut in concurrent.futures.as_completed(futures):
            repid, stats, pop = fut.result()
            statsDF = pd.DataFrame(stats, columns=DataRecord._fields)
            fixationsDF = pd.DataFrame(np.array(pop.fixations.array()))
            fixationsDF['repid'] = [repid]*len(fixationsDF.index)
            fixationsDF['ftime'] = pop.fixation_times
            with sqlite3.connect(args.outfile) as conn:
                statsDF.to_sql('data', conn, if_exists='append', index=False)
            with sqlite3.connect(args.fixationsfile) as conn:
                fixationsDF.to_sql(
                    'data', conn, if_exists='append', index=False)
