import fwdpy11
import fwdpy11.ts
import gzip
import sqlite3
import argparse
import sys
import os
import concurrent.futures
import libsequence.variant_matrix as vmatrix
import libsequence.summstats as sstats
import numpy as np
import pandas as pd
from collections import namedtuple

NLOCI = 10
LOCUS_BOUNDARIES = [(i, i + 11) for i in range(0, NLOCI * 11, 11)]

GenomeScanRecord = namedtuple(
    'GenomeScanRecord', ['generation', 'locus', 'window', 'pi', 'tajd', 'hprime'])

QtraitRecord = namedtuple('QtraitRecord', ['generation', 'zbar', 'wbar', 'vg'])


def make_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('suffix', type=str, help="Input file name suffix")
    parser.add_argument('lo', type=int, help="First replicate ID to process")
    parser.add_argument('hi', type=int, help="Last replicate ID to process")
    parser.add_argument('gscan_dbfile', type=str,
                        help="sqlite3 file name for genome scan output")
    parser.add_argument('qtrait_dbfile', type=str,
                        help="sqlite3 file name for quant get output")
    parser.add_argument('seed', type=int, default=None, help="RNG seed")

    parser.add_argument('--theta', type=float,
                        default=1e3, help="4Nu per locus")
    parser.add_argument('--nsam', type=int, default=50,
                        help="Sample size (number of diploids)")
    return parser


def process_replicate(argtuple):
    filename, repid, seed = argtuple

    with gzip.open(filename, 'rb') as f:
        pop = fwdpy11.SlocusPop.load_from_pickle_file(f)

    rng = fwdpy11.GSLrng(seed)

    # Add neutral mutations
    nm = fwdpy11.ts.infinite_sites(
        rng, pop, float(NLOCI) * args.theta / (4 * pop.N))

    # Get the node table for the pop and the
    # metadata for ancient samples
    nodes = np.array(pop.tables.nodes, copy=False)
    amd = np.array(pop.ancient_sample_metadata, copy=False)
    # These are the times for each ancient sample
    amt = nodes['time'][amd['nodes'][:, 0]]

    genome_scan = []
    qtraits = []
    for t in np.unique(amt):
        # Get the indexes of the metadata corresponding
        # to this time point
        sample_indexes_at_time = np.where(amt == t)[0]

        # Calc some simple stats about the overall pop'n
        mean_trait_value = amd['g'][sample_indexes_at_time].mean()
        vg = amd['g'][sample_indexes_at_time].var()
        wbar = amd['w'][sample_indexes_at_time].mean()
        qtraits.append(QtraitRecord(t, mean_trait_value, wbar, vg))

        # Get a random subset of the ancient samples
        # for our genome scan
        random_sample = np.random.choice(
            sample_indexes_at_time, args.nsam, replace=False)
        # Get the nodes associated w/those random samples
        samples = amd['nodes'][random_sample].flatten()

        # Simplify the tables w.r.to the samples.
        # Simplification results in a much smaller
        # set of trees to iterate over.
        tables, idmap = fwdpy11.ts.simplify(pop, samples)

        # These are the new indexes of our sample, after
        # simplification
        remapped = idmap[samples]

        # Now, we build a numpy array of the genotypes
        # containing both selected an neutral variants.
        # Because we use np.concatenate, to join
        # the neutral and selected sites, we have to be careful
        # to create windows where the sites are in the correct
        # order, which is increasing by position.  To
        # accomplish this, we rely on np.argsort
        gm = fwdpy11.ts.data_matrix_from_tables(
            tables, pop.mutations, remapped, True, True)
        pos = np.array(gm.neutral.positions + gm.selected.positions)
        sorted_pos_indexes = np.argsort(pos)
        all_sites = np.concatenate(
            (np.array(gm.neutral, copy=False), np.array(gm.selected, copy=False)))
        # We no longer need this object, and it is big.
        del gm
        # reorder the arrays
        pos = pos[sorted_pos_indexes]
        all_sites = all_sites[sorted_pos_indexes, :]

        # Now, we get a big libsequence object...
        vm = vmatrix.VariantMatrix(all_sites, pos)
        # ... and count the alleles
        ac = vm.count_alleles()

        # Delete object we no longer need
        del all_sites
        del vm
        # Now, we create our windows and iterate over them
        # via numpy "fancy" indexing methods:
        for locus, start_stop in enumerate(LOCUS_BOUNDARIES):
            for window, start in enumerate(range(*start_stop)):
                sites_in_window = np.where((pos >= start) & (
                    pos[sorted_pos_indexes] < start + 1.))[0]
                acw = ac[sites_in_window]
                pi = sstats.thetapi(acw)
                tajd = sstats.tajd(acw)
                hprime = sstats.hprime(acw, 0)
                genome_scan.append(GenomeScanRecord(
                    t, locus, window, pi, tajd, hprime))

    return repid, genome_scan, qtraits


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    np.random.seed(args.seed)
    used_seeds = []
    inputs = []
    for i in range(args.lo, args.hi + 1):
        fname = args.suffix + "/replicate{}.gz".format(i)
        if os.path.exists(fname) is False:
            raise ValueError("{} does not exist".format(fname))
        candidate_seed = np.random.randint(
            np.iinfo(np.uint32).max, dtype=np.uint32)
        while candidate_seed in used_seeds:
            candidate_seed = np.random.integer(
                np.iinfo(np.uint32).max, dtype=np.uint32)
        used_seeds.append(candidate_seed)
        inputs.append((fname, i, candidate_seed))

    used_seeds.clear()

    for i in [args.gscan_dbfile, args.qtrait_dbfile]:
        if os.path.exists(i):
            raise RuntimeError("File {} exists!".format(i))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_replicate, i) for i in inputs}
        for fut in concurrent.futures.as_completed(futures):
            repid, genome_scan, qtraits = fut.result()
            df = pd.DataFrame(genome_scan, columns=GenomeScanRecord._fields)
            df['repid'] = [repid] * len(df.index)
            with sqlite3.connect(args.gscan_dbfile) as conn:
                df.to_sql('data', conn, if_exists='append', index=False)
            genome_scan.clear()
            df = pd.DataFrame(qtraits, columns=QtraitRecord._fields)
            df['repid'] = [repid] * len(df.index)
            with sqlite3.connect(args.qtrait_dbfile) as conn:
                df.to_sql('data', conn, if_exists='append', index=False)
