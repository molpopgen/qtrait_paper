import fwdpy11
import fwdpy11.ts
import gzip
import sqlite3
import argparse
import sys
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
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

# NOTE: this implies that we'll be taking means of means later on!
LDRecord = namedtuple('LDRecord', ['generation', 'intralocus_zns', 'interlocus_zns',
                                   'intralocus_D', 'interlocus_D'])



def make_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('suffix', type=str, help="Input file name suffix")
    parser.add_argument('lo', type=int, help="First replicate ID to process")
    parser.add_argument('hi', type=int, help="Last replicate ID to process")
    parser.add_argument('gscan_dbfile', type=str,
                        help="sqlite3 file name for genome scan output")
    parser.add_argument('qtrait_dbfile', type=str,
                        help="sqlite3 file name for quant get output")
    parser.add_argument('ld_dbfile', type=str,
                        help="sqlite3 file name for quant get output")
    parser.add_argument('seed', type=int, default=None, help="RNG seed")

    parser.add_argument('--theta', type=float,
                        default=1e3, help="4Nu per locus")
    parser.add_argument('--nsam', type=int, default=50,
                        help="Sample size (number of diploids)")
    return parser


def pairwise_LD(t, pop, tables, samples, idmap):
    """
    Returns D and r^2 for all pairs of sites,
    along with site metadata
    """
    gm = fwdpy11.ts.data_matrix_from_tables(
        tables, pop.mutations, idmap[samples], False, True)
    selected = gm.selected
    nsites = selected.shape[0]
    nsam = selected.shape[1]

    pos = np.array(selected.positions)
    data = np.array(selected, copy=False)
    daf = np.sum(data, axis=1) / nsam
    daf_filter = np.where((daf >= 1e-2) & (daf < 1.0))[0]

    pos = pos[daf_filter]
    data = data[daf_filter, :]
    daf = daf[daf_filter]

    nsites = len(daf_filter)
    rsq_indexes = np.triu_indices(nsites, 1)

    rsq = np.ones(len(data[rsq_indexes]))
    D = np.zeros(len(data[rsq_indexes]))
    pos1 = np.zeros(len(rsq))
    pos2 = np.zeros(len(rsq))
    idx = 0
    for i in range(nsites - 1):
        for j in range(i + 1, nsites):
            p0 = daf[i]
            p1 = daf[j]
            temp = data[i, :] + data[j, :]
            p11 = len(np.where(temp == 2)[0]) / nsam
            pos1[idx] = pos[i]
            pos2[idx] = pos[j]
            D[idx] = p11 - p0 * p1
            rsq[idx] = np.power(D[idx], 2.0) / \
                (p0 * (1.0 - p0) * p1 * (1.0 - p1))
            idx += 1

    # assign a locus for all sites
    locus1 = np.zeros(len(rsq), dtype=np.int32)
    locus2 = np.zeros(len(rsq), dtype=np.int32)
    for locus, start_stop in enumerate(LOCUS_BOUNDARIES):
        for window, start in enumerate(range(*start_stop)):
            w1 = np.where((pos1 >= start) & (pos1 < start + 1.))[0]
            w2 = np.where((pos2 >= start) & (pos1 < start + 1.))[0]
            locus1[w1] = locus
            locus2[w2] = locus

    intralocus = np.zeros(len(rsq), dtype=np.int32)
    intralocus[np.where(locus1 == locus2)[0]] = 1
    intralocus_pairs = np.where(intralocus == 1)[0]
    interlocus_pairs = np.where(intralocus == 0)[0]
    rsq_intra = np.nan
    rsq_inter = np.nan
    D_intra = np.nan
    D_inter = np.nan
    if len(intralocus_pairs) > 0:
        rsq_intra = rsq[intralocus_pairs].mean()
        D_intra = D[intralocus_pairs].mean()

    if len(interlocus_pairs) > 0:
        rsq_inter = rsq[interlocus_pairs].mean()
        D_inter = D[interlocus_pairs].mean()

    return LDRecord(t, rsq_intra, rsq_inter, D_intra, D_inter)


def genome_scan_stats(t, mutations, tables, idmap, amd, sample_indexes_at_time):
    # Get a random subset of the ancient samples
    # for our genome scan
    random_sample = np.random.choice(
        sample_indexes_at_time, args.nsam, replace=False)
    # Get the nodes associated w/those random samples
    samples = amd['nodes'][random_sample].flatten()
    remapped = idmap[samples]

    # # Simplify the tables w.r.to the samples.
    # # Simplification results in a much smaller
    # # set of trees to iterate over.
    # tables, idmap = fwdpy11.ts.simplify(pop, samples)

    # # These are the new indexes of our sample, after
    # # simplification
    # remapped = idmap[samples]

    # Now, we build a numpy array of the genotypes
    # containing both selected an neutral variants.
    # Because we use np.concatenate, to join
    # the neutral and selected sites, we have to be careful
    # to create windows where the sites are in the correct
    # order, which is increasing by position.  To
    # accomplish this, we rely on np.argsort
    gm = fwdpy11.ts.data_matrix_from_tables(
        tables, mutations, remapped, True, True)

    pos = np.array(gm.neutral.positions + gm.selected.positions)
    sorted_pos_indexes = np.argsort(pos)
    all_sites = np.array(gm.neutral, copy=False)
    if len(gm.selected.positions) > 0:
        all_sites = np.concatenate((all_sites, np.array(gm.selected, copy=False)))
    # We no longer need this object, and it is big.
    # reorder the arrays
    pos = pos[sorted_pos_indexes]
    all_sites = all_sites[sorted_pos_indexes, :]

    # Now, we get a big libsequence object...
    vm = vmatrix.VariantMatrix(all_sites, pos)
    # ... and count the alleles
    ac = vm.count_alleles()
    # Now, we create our windows and iterate over them
    # via numpy "fancy" indexing methods:
    rv = []
    for locus, start_stop in enumerate(LOCUS_BOUNDARIES):
        for window, start in enumerate(range(*start_stop)):
            sites_in_window = np.where((pos >= start) & (
                pos[sorted_pos_indexes] < start + 1.))[0]
            acw = ac[sites_in_window]
            pi = sstats.thetapi(acw)
            tajd = sstats.tajd(acw)
            hprime = sstats.hprime(acw, 0)
            rv.append(GenomeScanRecord(
                t, locus, window, pi, tajd, hprime))

    return rv


def process_replicate(argtuple):
    filename, repid, seed = argtuple

    with gzip.open(filename, 'rb') as f:
        pop = fwdpy11.SlocusPop.load_from_pickle_file(f)

    rng = fwdpy11.GSLrng(seed)
    np.random.seed(seed)

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
    ld = []
    for t in np.unique(amt):
        # Get the indexes of the metadata corresponding
        # to this time point
        sample_indexes_at_time = np.where(amt == t)[0]

        # Calc some simple stats about the overall pop'n
        mean_trait_value = amd['g'][sample_indexes_at_time].mean()
        vg = amd['g'][sample_indexes_at_time].var()
        wbar = amd['w'][sample_indexes_at_time].mean()
        qtraits.append(QtraitRecord(t, mean_trait_value, wbar, vg))

        samples = amd['nodes'][sample_indexes_at_time].flatten()
        tables, idmap = fwdpy11.ts.simplify(pop, samples)

        ld.append(pairwise_LD(t, pop, tables, samples, idmap))
        genome_scan.extend(genome_scan_stats(
            t, pop.mutations, tables, idmap, amd, sample_indexes_at_time))

    return repid, genome_scan, qtraits, ld


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

    for i in [args.gscan_dbfile, args.qtrait_dbfile, args.ld_dbfile]:
        if os.path.exists(i):
            raise RuntimeError("File {} exists!".format(i))

    with concurrent.futures.ProcessPoolExecutor(max_workers=64) as executor:
        futures = {executor.submit(process_replicate, i) for i in inputs}
        for fut in concurrent.futures.as_completed(futures):
            repid, genome_scan, qtraits, ld = fut.result()
            df = pd.DataFrame(genome_scan, columns=GenomeScanRecord._fields)
            df['repid'] = [repid] * len(df.index)
            with sqlite3.connect(args.gscan_dbfile) as conn:
                df.to_sql('data', conn, if_exists='append', index=False)
            genome_scan.clear()
            df = pd.DataFrame(qtraits, columns=QtraitRecord._fields)
            df['repid'] = [repid] * len(df.index)
            with sqlite3.connect(args.qtrait_dbfile) as conn:
                df.to_sql('data', conn, if_exists='append', index=False)
            df = pd.DataFrame(ld, columns=LDRecord._fields)
            df['repid'] = [repid] * len(df.index)
            with sqlite3.connect(args.ld_dbfile) as conn:
                df.to_sql('data', conn, if_exists='append', index=False)

