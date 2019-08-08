import fwdpy11
import gzip
import sys
import argparse
import numpy as np
import pandas as pd
import sqlite3
from collections import namedtuple

Datum = namedtuple('Datum', ['generation', 'position', 'origin', 'ftime',
                             'esize', 'nsamples', 'W', 'G',
                             'Wbar', 'Gbar', 'VG'])


def make_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required arguments")
    required.add_argument("--infile", type=str, default=None,
                          help="Input file name")
    required.add_argument("--outfile", type=str, default=None,
                          help="Output file name")
    return parser


def getpop(args):
    with gzip.open(args.infile, 'rb') as f:
        pop = fwdpy11.DiploidPopulation.load_from_pickle_file(f)
    return pop


def process_pop(pop):
    data = []
    nodes = np.array(pop.tables.nodes, copy=False)
    # Get the metadata for ancient samples
    amd = np.array(pop.ancient_sample_metadata, copy=False)
    # Get the ancient sample time points
    times = nodes['time'][amd['nodes'][:, 0]]
    utimes = np.unique(times)

    fpos = []
    ftimes = []
    origins = []
    sites = np.array(pop.tables.sites, copy=False)
    muts = np.array(pop.tables.mutations, copy=False)
    for i, j in zip(pop.fixation_times, pop.fixations):
        if i > 10*pop.N and j.g <= 10*pop.N + 200:
            fpos.append(j.pos)
            ftimes.append((i, j.pos))
            origins.append((j.g, j.pos))
    ftimes = np.array([i[0] for i in sorted(ftimes, key=lambda x: x[1])])
    origins = np.array([i[0] for i in sorted(origins, key=lambda x: x[1])])
    fpos = np.array(sorted(fpos))
    for ut in utimes:
        mdidx = np.where(times == ut)[0]
        samples = amd['nodes'][mdidx].flatten()
        tables, idmap = fwdpy11.simplify_tables(pop.tables, samples)
        sites = np.array(tables.sites, copy=False)
        muts = np.array(tables.mutations, copy=False)
        tv = fwdpy11.TreeIterator(tables, idmap[samples], update_samples=True)
        for t in tv:
            l = t.left
            r = t.right
            f = np.where((fpos >= l) & (fpos < r))[0]
            for i in f:
                idx = np.where(sites['position'] == fpos[i])[0]
                if len(idx) > 0:
                    mut_idx = np.where(muts['site'] == idx)[0]
                    assert len(mut_idx == 1), "Bad mutation table error"
                    sbelow = t.samples_below(muts['node'][mut_idx])
                    individuals = np.unique(sbelow // 2)
                    mg = amd['g'][mdidx[individuals]].mean()
                    gbar = amd['g'][mdidx].mean()
                    vg = amd['g'][mdidx].var()
                    mw = amd['w'][mdidx[individuals]].mean()
                    wbar = amd['w'][mdidx].mean()
                    esize = pop.mutations[muts['key'][mut_idx[0]]].s
                    data.append(Datum(ut, fpos[i], origins[i],
                                      ftimes[i], esize,
                                      len(sbelow), mw, mg, wbar, gbar, vg))
    return data


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    if args.outfile is None:
        raise ValueError("output file name cannot be None")
    pop = getpop(args)
    data = process_pop(pop)
    df = pd.DataFrame(data, columns=Datum._fields)
    with sqlite3.connect(args.outfile) as conn:
        df.to_sql('data', conn, if_exists="replace")
