import fwdpy11
import gzip
import sys
import argparse
import numpy as np
import pandas as pd
import sqlite3
import libsequence
from collections import namedtuple

SweepData = namedtuple('SweepData', ['generation', 'locus', 'D', 'Hp', 'hdiv'])

LOCUS_BOUNDARIES = [(i, i + 11) for i in range(0, 10 * 11, 11)]
SREGIONS = [(i[0] + 5, i[0] + 6) for i in LOCUS_BOUNDARIES]
THETA_PER_LOCUS = 1e3
THETA = THETA_PER_LOCUS * len(LOCUS_BOUNDARIES)
NSAM = 50


def make_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required arguments")
    required.add_argument("--infile", type=str, default=None,
                          help="Input file name")
    required.add_argument("--outfile", type=str, default=None,
                          help="Output file name")
    required.add_argument('--seed', type=int,
                          default=None, help="RNG seed")
    return parser


def getpop(args):
    with gzip.open(args.infile, 'rb') as f:
        pop = fwdpy11.DiploidPopulation.load_from_pickle_file(f)
    return pop


def merge_matrix(dm):
    if dm.selected.shape[0] == 0:
        return np.array(dm.neutral_positions), np.array(dm.neutral)
    combined = np.concatenate((dm.neutral, dm.selected), axis=0)
    # combined the positions from the 2 inputs
    p = np.concatenate((dm.neutral_positions, dm.selected_positions))
    pindex = [i for i in range(len(p))]
    pindex = sorted(pindex, key=lambda x: p[x])
    # We return the sorted positions
    # and the columns re-sorted according to position
    return p[pindex], combined[pindex, :]


def process_pop(pop):
    nodes = np.array(pop.tables.nodes, copy=False)
    # Get the metadata for ancient samples
    amd = np.array(pop.ancient_sample_metadata, copy=False)
    # Get the ancient sample time points
    times = nodes['time'][amd['nodes'][:, 0]]
    utimes = np.unique(times)

    sweep_data = []
    rng = fwdpy11.GSLrng(args.seed)
    np.random.seed(args.seed)
    fwdpy11.infinite_sites(rng, pop, THETA/(4*pop.N))
    for ut in utimes:
        mdidx = np.where(times == ut)[0]
        mdidx = np.random.choice(mdidx, NSAM, replace=False)
        samples = amd['nodes'][mdidx].flatten()
        tables, idmap = fwdpy11.simplify_tables(pop.tables, samples)
        for locus, dm in enumerate(fwdpy11.DataMatrixIterator(tables,
                                                              idmap[samples],
                                                              SREGIONS, True, True)):
            pos, mat = merge_matrix(dm)
            vm = libsequence.VariantMatrix(mat, pos)
            ac = vm.count_alleles()
            D = libsequence.tajd(ac)
            Hp = libsequence.hprime(ac, 0)
            hdiv = libsequence.haplotype_diversity(vm)
            sweep_data.append(SweepData(int(ut), locus, D, Hp, hdiv))

    return sweep_data


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    if args.outfile is None:
        raise ValueError("output file name cannot be None")
    pop = getpop(args)
    sweep_data = process_pop(pop)
    df_sweep = pd.DataFrame(sweep_data, columns=SweepData._fields)
    with sqlite3.connect(args.outfile) as conn:
        df_sweep.to_sql('data_sweep', conn, if_exists='replace')
