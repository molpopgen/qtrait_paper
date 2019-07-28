import fwdpy11
import argparse
import sys
import gzip
import numpy as np


def make_parser():
    """
    Create a command-line interface to the script.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required arguments")
    required.add_argument("--popsize", "-N", type=int,
                          help="Diploid population size")
    required.add_argument("--filename", '-f', type=str, default=None,
                          help="Output file name.")
    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument("--seed", type=int, default=42,
                          help="Random number seed")
    optional.add_argument("--mu", "-m", type=float, default=1e-3,
                          help="Mutation rate (per gamete, per generation)")
    optional.add_argument("--sigma", "-s", type=float, default=0.25,
                          help="Standard deviation of Gaussian"
                          "distribution of mutational effects")
    optional.add_argument("--rho", type=float,
                          default=1e3, help="4Nr per locus")
    optional.add_argument("--opt", type=float, default=1.0,
                          help="Value of new phenotypic optimum")
    optional.add_argument('--nloci', type=int, default=10,
                          help="Number of independently-assorting loci")

    return parser


class Recorder(object):
    def __init__(self, N):
        self.individuals = np.arange(N, dtype=np.int32)
        self.N = N

    def __call__(self, pop, recorder):
        if pop.generation >= 10 * self.N:
            if pop.generation <= 10*self.N + 200:
                recorder.assign(self.individuals)


def runsim(args):
    popsizes = np.array([args.popsize]*20*args.popsize, dtype=np.int32)
    locus_boundaries = [(i, i + 11) for i in range(0, args.nloci * 11, 11)]
    pop = fwdpy11.DiploidPopulation(args.popsize, locus_boundaries[-1][1])
    sregions = [fwdpy11.GaussianS(i[0] + 5, i[0] + 6, 1, args.sigma)
                for i in locus_boundaries]
    recregions = [fwdpy11.PoissonInterval(
        *i, args.rho / (4 * args.popsize)) for i in locus_boundaries]
    recregions.extend([fwdpy11.BinomialPoint(i[1], 0.5)
                       for i in locus_boundaries[:-1]])

    optima = fwdpy11.GSSmo([(0, 0, 1), (10*args.popsize, args.opt, 1)])
    p = {'nregions': [],  # No neutral mutations -- add them later!
         'gvalue': fwdpy11.Additive(2.0, optima),
         'sregions': sregions,
         'recregions': recregions,
         'rates': (0.0, args.mu, None),
         # Keep mutations at frequency 1 in the pop if they affect fitness.
         'prune_selected': False,
         'demography':  np.array(popsizes, dtype=np.uint32)
         }
    params = fwdpy11.ModelParams(**p)
    rng = fwdpy11.GSLrng(args.seed)
    s = Recorder(args.popsize)
    fwdpy11.evolvets(rng, pop, params, 100, s,
                     suppress_table_indexing=True,
                     track_mutation_counts=True)
    return pop


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    if args.filename is None:
        raise ValueError("output file name cannot be None")
    pop = runsim(args)
    with gzip.open(args.filename, 'wb') as f:
        pop.pickle_to_file(f)
