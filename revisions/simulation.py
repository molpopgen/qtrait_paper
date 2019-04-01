"""
Requires fwdpy11 0.3.0.

This is a "hack" of the simulation found in fwdpy11/examples/gss

NOTE: the usage is a bit wacko.  I allow concurrent futures to
automate many replicates or running one replicate with a known 
replicate ID.  These options make the CLI messy.

Evolution of a single genomic region under Gaussian Stabilizing Selection
with a moving optimum, or "GSSmo".

The simulation allows the recording of quantitative-genetics statistics
over time as well as the preservation of samples into the tree sequence
after the optimum shift.

The demographic model is constant-size Wright-Fisher for 10N generations
around an optimal trait value of zero.  Then, the optimum changes and
the population continues to evolve for a length of time specified
by the user.

The mutation model is uniform across the genome with effect sizes
given by a Gaussian distribution with mean zero and a user-specified
standard deviation.

The genetic model is additive effects with fitness determined by
Gaussian stabilizing selection based on the squared distance from
the optimum.
"""

import fwdpy11 as fp11
import fwdpy11.model_params
import fwdpy11.genetic_values
import fwdpy11.ts
import fwdpy11.tsrecorders
import fwdpy11.wright_fisher_ts
import math
import numpy as np
import gzip
import sys
import argparse
import concurrent.futures
import pandas as pd
import sqlite3
from collections import defaultdict
from collections import namedtuple

LDRecord = namedtuple('LDRecord', ['generation', 'intralocus_zns', 'interlocus_zns',
                                   'intralocus_D', 'interlocus_D'])


def make_parser():
    """
    Create a command-line interface to the script.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required arguments")
    required.add_argument("--popsize", "-N", type=int,
                          help="Diploid population size")
    required.add_argument("--mu", "-m", type=float,
                          help="Mutation rate (per gamete, per generation)")
    required.add_argument("--filename", "-f", type=str, default=None,
                          help="Prefix of output file name.  The population will be"
                          "pickled to this file, and the file will be"
                          "compressed using the gzip algorithm")
    required.add_argument("--ldfile", type=str, default=None,
                          help="sqlite3 database for LD data")

    optional = parser.add_argument_group("Optional arguments")
    required.add_argument("--sigma", "-s", type=float, default=None,
                          help="Standard deviation of Gaussian"
                          "distribution of mutational effects")
    optional.add_argument("--rho", type=float, default=1000.0,
                          help="Scaled recombination rate, rho=4Nr")
    optional.add_argument("--VS", type=float, default=1.0,
                          help="Inverse strength of stabilizing selection")
    optional.add_argument("--opt", type=float, default=1.0,
                          help="Value of new phenotypic optimum")
    optional.add_argument('--nloci', type=int, default=10,
                          help="Number of independently-assorting loci")
    optional.add_argument("--time", type=float, default=0.1,
                          help="Amount of time to simulate past"
                          "optimum shift, in units of N")
    optional.add_argument("--preserve", type=int, default=-1,
                          help="Record ancient samples every X generations"
                          "after optimum shift. A value of -1 means not"
                          "to record.")
    optional.add_argument("--num_ind", '-n', type=int, default=50,
                          help="Number of diploids to record as"
                          "ancient samples")
    optional.add_argument("--seed", type=int, default=42,
                          help="Random number seed.")
    optional.add_argument("--repid", type=int,
                          default=None, help="Replicate ID")
    optional.add_argument("--nreps", type=int, default=1,
                          help="Number of replicates to run.")
    optional.add_argument("--max_workers", type=int, default=None,
                          help="Max number of processes to run at once.")

    return parser


def validate_arguments(args):
    """
    Validate input arguments.
    Note: this is likely incomplete.
    """
    if args.popsize < 0:
        raise ValueError("Population size must be non-negative")
    if args.mu < 0 or math.isfinite(args.mu) is False:
        raise ValueError("Mutation rate must be non-negative and finite")
    if args.sigma is not None and (args.sigma < 0 or math.isfinite(args.sigma) is False):
        raise ValueError("Std. dev. of distribution of effect sizes"
                         "must be non-negative and finite")
    if args.filename is None:
        raise ValueError("ouptut filename cannot be None")

    if args.preserve > 0:
        if args.num_ind > args.popsize:
            raise ValueError("Number of ancient samples is"
                             "greater than population size")


class Recorder(object):
    """
    fwdpy11 allows you to define objects that record data
    from populations during simulation.  Such objects must
    be callable, and the easiest way to do things is to
    create a class with a __call__ function.
    """

    def __init__(self, maxtime, interval, nsam, nloci):
        self.data = []
        self.maxtime = maxtime
        self.interval = interval
        self.nsam = nsam
        self.locus_boundaries = [(i, i + 11)
                                 for i in range(0, args.nloci * 11, 11)]
        self.ld = []

    def getld(self, pop):
        # Get the LD
        keys = defaultdict(lambda: 0)
        for g in pop.gametes:
            if g.n > 0:
                for k in g.smutations:
                    keys[k] += g.n

        segregating = [key for key, value in keys.items() if value >
                       0 and value < 2 * pop.N]
        pos = np.array(sorted([pop.mutations[i].pos for i in segregating]))
        sorted_key_map = {}
        for i in segregating:
            idx = np.where(pos == pop.mutations[i].pos)[0][0]
            sorted_key_map[i] = idx

        rsq_indexes = np.triu_indices(len(pos), 1)
        rsq = np.ones(len(rsq_indexes[0]))
        D = np.ones(len(rsq_indexes[0]))
        rsq.fill(np.nan)
        D.fill(np.nan)

        genotypes = np.zeros(len(pos) * 2 * pop.N).reshape(len(pos), 2 * pop.N)

        for i, d in enumerate(pop.diploids):
            for k in pop.gametes[d.first].smutations:
                if k in sorted_key_map:
                    genotypes[sorted_key_map[k], 2 * i] = 1
            for k in pop.gametes[d.second].smutations:
                if k in sorted_key_map:
                    genotypes[sorted_key_map[k], 2 * i + 1] = 1

        daf = np.sum(genotypes, axis=1) / genotypes.shape[1]
        pos1 = np.zeros(len(rsq))
        pos2 = np.zeros(len(rsq))
        pos1.fill(np.nan)
        pos2.fill(np.nan)
        idx = 0
        for i in range(genotypes.shape[0] - 1):
            for j in range(i + 1, genotypes.shape[0]):
                p0 = daf[i]
                p1 = daf[j]
                pos1[idx] = pos[i]
                pos2[idx] = pos[j]
                if p0 >= 1e-2 and p1 >= 1e-2:
                    temp = genotypes[i, :] + genotypes[j, :]
                    p11 = len(np.where(temp == 2)[0]) / genotypes.shape[1]
                    D[idx] = p11 - p0 * p1
                    rsq[idx] = np.power(D[idx], 2.0) / \
                        (p0 * (1.0 - p0) * p1 * (1.0 - p1))
                idx += 1

        locus1 = np.zeros(len(rsq), dtype=np.int32)
        locus2 = np.zeros(len(rsq), dtype=np.int32)
        for locus, start_stop in enumerate(self.locus_boundaries):
            for window, start in enumerate(range(*start_stop)):
                w1 = ((pos1 >= start) & (pos1 < start + 1.)).nonzero()
                w2 = ((pos2 >= start) & (pos2 < start + 1.)).nonzero()
                locus1[w1] = locus
                locus2[w2] = locus

        intralocus = np.zeros(len(rsq), dtype=np.int32)
        intralocus[np.where(locus1 == locus2)[0]] = 1
        intralocus_pairs = np.where(intralocus == 1)[0]
        interlocus_pairs = np.where(intralocus == 0)[0]
        if len(intralocus_pairs) > 0:
            rsq_intra = np.nanmean(rsq[intralocus_pairs])
            D_intra = np.nanmean(D[intralocus_pairs])
        else:
            rsq_intra = np.nan
            D_intra = np.nan

        if len(interlocus_pairs) > 0:
            rsq_inter = np.nanmean(rsq[interlocus_pairs])
            D_inter = np.nanmean(D[interlocus_pairs])
        else:
            rsq_inter = np.nan
            D_inter = np.nan

        return rsq_intra, rsq_inter, D_intra, D_inter

    def __call__(self, pop, recorder):
        if self.interval > 0 and pop.generation >= 10 * pop.N:
            if pop.generation <= 10 * pop.N + self.maxtime and pop.generation % self.interval == 0.0:
                if self.nsam < pop.N:
                    s = np.random.choice(pop.N, self.nsam, replace=False)
                else:
                    s = np.arange(pop.N, dtype=np.int32)

                recorder.assign(s)
        if pop.generation >= 10 * pop.N and pop.generation <= 10 * pop.N + 4 * pop.N:
            ld = self.getld(pop)
            self.ld.append(LDRecord(pop.generation, *ld))


def runsim(argtuple):
    """
    Run the simulation and deliver output to files.
    """
    args, repid, seed = argtuple
    locus_boundaries = [(i, i + 11) for i in range(0, args.nloci * 11, 11)]
    pop = fp11.SlocusPop(args.popsize, locus_boundaries[-1][1])
    sregions = [fwdpy11.GaussianS(i[0] + 5, i[0] + 6, 1, args.sigma)
                for i in locus_boundaries]
    recregions = [fwdpy11.PoissonInterval(
        *i, args.rho / (4 * args.popsize)) for i in locus_boundaries]
    recregions.extend([fwdpy11.BinomialPoint(i[1], 0.5)
                       for i in locus_boundaries[:-1]])

    rng = fp11.GSLrng(seed)

    GSS = fp11.genetic_values.GSS(opt=0, VS=1)

    popsizes = [args.popsize] * (10 * args.popsize - 1)
    p = {'nregions': [],  # No neutral mutations -- add them later!
         'gvalue': fwdpy11.genetic_values.SlocusAdditive(2.0, GSS),
         'sregions': sregions,
         'recregions': recregions,
         'rates': (0.0, args.mu, None),
         # Keep mutations at frequency 1 in the pop if they affect fitness.
         'prune_selected': False,
         'demography':  np.array(popsizes, dtype=np.uint32)
         }
    params = fp11.model_params.ModelParams(**p)

    r = Recorder(int(args.time * float(args.popsize)),
                 args.preserve, args.num_ind, args.nloci)
    fwdpy11.wright_fisher_ts.evolve(
        rng, pop, params, 100, r, suppress_table_indexing=True,
        remove_extinct_variants=False)

    GSSmo = fp11.genetic_values.GSSmo(
        [(pop.generation, 0, args.VS), (10 * args.popsize, args.opt, args.VS)])
    p['gvalue'] = fwdpy11.genetic_values.SlocusAdditive(2.0, GSSmo)
    p['demography'] = np.array(
        [args.popsize] * (1 + 10 * args.popsize), dtype=np.uint32)
    params = fp11.model_params.ModelParams(**p)
    fwdpy11.wright_fisher_ts.evolve(
        rng, pop, params, 100, r, suppress_table_indexing=True,
        track_mutation_counts=True)

    fname = args.filename + "{}.gz".format(repid)
    with gzip.open(fname, 'wb') as f:
        pop.pickle_to_file(f)

    return repid, r.ld


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    validate_arguments(args)

    if args.repid is not None:
        print("running with fixed replicate id")
        repid, ld = runsim((args, args.repid, args.seed))
        df = pd.DataFrame(ld, columns = LDRecord._fields)
        df['repid'] = [repid]*len(df.index)
        with sqlite3.connect(args.ldfile) as conn:
            df.to_sql('data', conn, index=False)
    else:
        np.random.seed(args.seed)

        seeds = []
        for i in range(args.nreps):
            candidate = np.random.randint(
                np.iinfo(np.uint32).max, dtype=np.uint32)
            while candidate in seeds:
                candidate = np.random.randint(
                    np.iinfo(np.uint32).max, dtype=np.uint32)
            seeds.append(candidate)

        with concurrent.futures.ProcessPoolExecutor(max_workers=args.max_workers) as executor:
            futures = {executor.submit(
                runsim, (args, i[0], i[1])): i for i in enumerate(seeds)}

            for fut in concurrent.futures.as_completed(futures):
                repid, ld = fut.result()
                df = pd.DataFrame(ld, columns = LDRecord._fields)
                df['repid'] = [repid]*len(df.index)
                with sqlite3.connect(args.ldfile) as conn:
                    df.to_sql('data', conn, index=False)

