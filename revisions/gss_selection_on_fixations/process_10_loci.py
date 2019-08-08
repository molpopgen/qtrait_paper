import fwdpy11
import libsequence
import gzip
import sys
import sqlite3
import numpy as np
import pandas as pd
from collections import namedtuple

NLOCI = 10
LOCUS_LENGTH = 11
THETA_PER_LOCUS = 1e3
DIPLOID_SAMPLESIZE = 50

STEPSIZE = 0.1
WINDOW_LEFTS = np.arange(0, NLOCI*LOCUS_LENGTH, STEPSIZE)
WINDOWS = [(i, i+STEPSIZE) for i in WINDOW_LEFTS]

Datum = namedtuple('Datum', ['left_edge', 'locus', 'window', 'pwindow',
                             'pi', 'D', 'Hp', 'nhaps', 'hdiv'])


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


def get_stats(tables, samples):
    tables_for_sample, idmap = fwdpy11.simplify_tables(tables, samples)
    stats = []
    for k, l in zip(enumerate(fwdpy11.DataMatrixIterator(tables_for_sample,
                                                         idmap[samples],
                                                         WINDOWS, True, True)),
                    WINDOW_LEFTS):
        i, dm = k
        locus = int(i*STEPSIZE / LOCUS_LENGTH)
        window = i % (LOCUS_LENGTH/STEPSIZE)
        window_in_paper = int(window*STEPSIZE)
        pos, data = merge_matrix(dm)
        vm = libsequence.VariantMatrix(data, pos)
        ac = vm.count_alleles()
        pi = libsequence.thetapi(ac)
        D = libsequence.tajd(ac)
        Hp = libsequence.hprime(ac, 0)
        nhaps = libsequence.number_of_haplotypes(vm)
        hdiv = libsequence.haplotype_diversity(vm)
        stats.append(Datum(l, locus, window, window_in_paper, pi, D, Hp, nhaps, hdiv))

    return stats


def process_file(rng, filename):
    with gzip.open(filename, 'rb') as f:
        pop = fwdpy11.DiploidPopulation.load_from_pickle_file(f)

    fwdpy11.infinite_sites(rng, pop, THETA_PER_LOCUS*NLOCI/(4*pop.N))
    nodes = np.array(pop.tables.nodes, copy=False)
    # Get the metadata for ancient samples
    amd = np.array(pop.ancient_sample_metadata, copy=False)
    # Get the ancient sample time points
    times = nodes['time'][amd['nodes'][:, 0]]
    utimes = np.unique(times)

    thresholds = [0.1, 0.5, 0.9]
    dflist = []
    for ut in utimes:
        mdidx = np.where(times == ut)[0]
        gbar = amd['g'][mdidx].mean()
        if gbar >= thresholds[0]:
            individuals = np.random.choice(
                mdidx, DIPLOID_SAMPLESIZE, replace=False)
            sample_nodes = amd['nodes'][individuals].flatten()
            stats = get_stats(pop.tables, sample_nodes)
            df = pd.DataFrame(stats, columns=Datum._fields)
            df['generation'] = [ut]*len(df.index)
            df['gbar'] = [gbar]*len(df.index)
            df['threshold'] = [thresholds[0]]*len(df.index)
            thresholds.pop(0)
            dflist.append(df)
            if len(thresholds) == 0:
                break

    df = pd.concat(dflist)
    return pop, df


def fixations(pop):
    ftimes = np.array(pop.fixation_times, dtype=np.int32)
    esizes = np.zeros(len(ftimes))
    otimes = np.zeros(len(ftimes), dtype=np.int32)
    position = np.zeros(len(ftimes))
    for i, j in enumerate(pop.fixations):
        esizes[i] = j.s
        otimes[i] = j.g
        position[i] = j.pos
    return pd.DataFrame({'position': position, 'origin': otimes,
                         'fixation': ftimes, 'esize': esizes})


if __name__ == "__main__":
    filename = sys.argv[1]
    dbname = sys.argv[2]
    seed = int(sys.argv[3])
    np.random.seed(seed)
    rng = fwdpy11.GSLrng(seed)
    pop, results = process_file(rng, filename)
    fixation_df = fixations(pop)

    with sqlite3.connect(dbname) as conn:
        results.to_sql('data', conn, if_exists='replace')
        fixation_df.to_sql('fixations', conn, if_exists='replace')
