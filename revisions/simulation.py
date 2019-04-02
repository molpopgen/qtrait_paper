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

# LDRecord = namedtuple('LDRecord', ['generation', 'intralocus_zns', 'interlocus_zns',
#                                    'intralocus_D', 'interlocus_D'])
LDRecord = namedtuple(
    'LDRecord', ['generation', 'pos1', 'pos2', 'c1', 'c2', 'e1', 'e2', 'D', 'rsq'])


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

    def __init__(self, maxtime, interval, nsam, dbname, repid):
        self.data = []
        self.maxtime = maxtime
        self.interval = interval
        self.nsam = nsam
        self.ld = []
        self.dbname = dbname
        self.repid = repid

    def getld(self, pop):
        mc = np.array(pop.mcounts)
        segregating = ((mc > 0) & (mc < 2 * pop.N)).nonzero()[0]
        sorted_key_map = defaultdict(lambda: 0)
        for i, j in enumerate(segregating):
            sorted_key_map[j] = i
        pos = np.array([pop.mutations[i].pos for i in segregating])
        genotypes = np.zeros(len(pos) * 2 * pop.N).reshape(len(pos), 2 * pop.N)

        for i, d in enumerate(pop.diploids):
            for k in pop.gametes[d.first].smutations:
                if k in sorted_key_map:
                    genotypes[sorted_key_map[k], 2 * i] = 1
            for k in pop.gametes[d.second].smutations:
                if k in sorted_key_map:
                    genotypes[sorted_key_map[k], 2 * i + 1] = 1

        dcounts = np.sum(genotypes, axis=1).astype(np.int32)
        idx = 0
        for i in range(genotypes.shape[0] - 1):
            for j in range(i + 1, genotypes.shape[0]):
                p0 = dcounts[i] / (2 * pop.N)
                p1 = dcounts[j] / (2 * pop.N)
                temp = genotypes[i, :] + genotypes[j, :]
                p11 = len(np.where(temp == 2)[0]) / genotypes.shape[1]
                D = p11 - p0 * p1
                rsq = np.power(D, 2.0) / \
                    (p0 * (1.0 - p0) * p1 * (1.0 - p1))
                idx += 1
                self.ld.append(LDRecord(pop.generation,
                                        pos[i], pos[j], dcounts[i], dcounts[j],
                                        pop.mutations[segregating[i]].s,
                                        pop.mutations[segregating[j]].s,
                                        D, rsq))
        if len(self.ld) > 200000:
            temp = pd.DataFrame(self.ld, columns=LDRecord._fields)
            temp['repid'] = np.array(
                [self.repid] * len(temp.index), dtype=np.int32)
            with sqlite3.connect(self.dbname) as conn:
                temp.to_sql('data', conn, if_exists='append',
                            chunksize=50000, index=False)
            self.ld.clear()

    def __call__(self, pop, recorder):
        if self.interval > 0 and pop.generation >= 10 * pop.N:
            if pop.generation <= 10 * pop.N + self.maxtime and pop.generation % self.interval == 0.0:
                if self.nsam < pop.N:
                    s = np.random.choice(pop.N, self.nsam, replace=False)
                else:
                    s = np.arange(pop.N, dtype=np.int32)

                recorder.assign(s)
        if pop.generation >= 10 * pop.N and pop.generation <= 10 * pop.N + 4 * pop.N:
            self.getld(pop)


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

    dbname = args.filename + "{}_ld.sqlite3".format(repid)
    r = Recorder(int(args.time * float(args.popsize)),
                 args.preserve, args.num_ind, dbname, repid)
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

    if len(r.ld) > 0:
        df = pd.DataFrame(r.ld, columns=LDRecord._fields)
        with sqlite3.connect(dbname) as conn:
            df.to_sql('data', conn, if_exists='append', index=False,
                      chunksize=50000)
        r.ld.clear()

    fname = args.filename + "{}.gz".format(repid)
    with gzip.open(fname, 'wb') as f:
        pop.pickle_to_file(f)

    # df = pd.DataFrame(r.ld, columns=LDRecord._fields)
    # df['repid'] = [repid] * len(df.index)

    # with sqlite3.connect(fname) as conn:
    #     df.to_sql('data', conn, if_exists='append', index=False)

    # return repid, r.ld


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    validate_arguments(args)

    if args.repid is not None:
        print("running with fixed replicate id")
        rv = runsim((args, args.repid, args.seed))
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
                rv = fut.result()
