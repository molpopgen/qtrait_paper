"""
Requires fwdpy11 0.3.0.

This is a "hack" of the simulation found in fwdpy11/examples/gss

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

    def __init__(self, maxtime, interval, nsam):
        self.data = []
        self.maxtime = maxtime
        self.interval = interval
        self.nsam = nsam

    def __call__(self, pop, recorder):
        if self.interval > 0 and pop.generation >= 10 * pop.N:
            if pop.generation <= 10 * pop.N + self.maxtime and pop.generation % self.interval == 0.0:
                if self.nsam < pop.N:
                    s = np.random.choice(pop.N, self.nsam, replace=False)
                else:
                    s = np.arange(pop.N, dtype=np.int32)

                recorder.assign(s)


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
                 args.preserve, args.num_ind)
    fwdpy11.wright_fisher_ts.evolve(
        rng, pop, params, 100, r, suppress_table_indexing=True)

    GSSmo = fp11.genetic_values.GSSmo(
        [(pop.generation, 0, args.VS), (10 * args.popsize, args.opt, args.VS)])
    p['gvalue'] = fwdpy11.genetic_values.SlocusAdditive(2.0, GSSmo)
    p['demography'] = np.array(
        [args.popsize] * (1 + 10 * args.popsize), dtype=np.uint32)
    params = fp11.model_params.ModelParams(**p)
    fwdpy11.wright_fisher_ts.evolve(
        rng, pop, params, 100, r, suppress_table_indexing=False,
        track_mutation_counts=True)

    fname = args.filename + "{}.gz".format(repid)
    with gzip.open(fname, 'wb') as f:
        pop.pickle_to_file(f)

    return True


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    validate_arguments(args)

    np.random.seed(args.seed)

    seeds = []
    for i in range(args.nreps):
        candidate = np.random.randint(np.iinfo(np.uint32).max, dtype=np.uint32)
        while candidate in seeds:
            candidate = np.random.randint(0, dtype=np.uint32)
        seeds.append(candidate)

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {executor.submit(
            runsim, (args, i[0], i[1])): i for i in enumerate(seeds)}

        for fut in concurrent.futures.as_completed(futures):
            result = fut.result()
            if result is not True:
                raise RuntimeError("Unexpected return value from process")
