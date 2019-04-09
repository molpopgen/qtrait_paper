import sqlite3
import fwdpy11
import glob
import argparse
import sys
import os
import gzip
import pandas as pd
import concurrent.futures
from collections import namedtuple

Datum = namedtuple('Datum',['s','g','ftime'])

def make_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('stub',type=str)
    parser.add_argument('lo',type=int)
    parser.add_argument('hi',type=int)
    parser.add_argument('outfile',type=str)
    parser.add_argument('cores',type=int)

    return parser

def process_file(args):
    fname, repid, mu = args
    with gzip.open(fname, 'rb') as f:
        pop = fwdpy11.SlocusPop.load_from_pickle_file(f)

    fixations = []
    for i,j in zip(pop.fixations, pop.fixation_times):
        fixations.append(Datum(i.s, i.g, j))

    return fixations, repid, mu


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    if os.path.exists(args.outfile):
        raise RuntimeError("{} already exists".format(args.outfile))

    model_params = {'lomu': 2.5e-4,
            'midmu': 1e-3,
            'himu': 5e-3}

    inputs = []
    for i in range(args.lo, args.hi+1):
        for key, value in model_params.items():
            fname=args.stub + "/{}/replicate{}.gz".format(key, i)
            if os.path.exists(fname) is False:
                raise ValueError("{} does not exist".format(fname))
            inputs.append((fname, i, value))

    with concurrent.futures.ProcessPoolExecutor(max_workers=int(args.cores/2)) as executor:
        futures = {executor.submit(process_file, i) for i in inputs}

        for fut in concurrent.futures.as_completed(futures):
            fixations, repid, mu = fut.result()
            df = pd.DataFrame(fixations, columns=Datum._fields)
            df['repid'] = [repid]*len(df.index)
            df['mu'] = [mu]*len(df.index)

            with sqlite3.connect(args.outfile) as conn:
                df.to_sql('data', conn, if_exists='append', index=False)
