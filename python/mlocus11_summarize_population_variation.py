import sys
import os
import pickle
import tarfile
import argparse
import lzma
import concurrent.futures
import numpy as np
import pandas as pd
import sqlite3
import os
from fwdpy11.sampling import sample_separate
from fwdpy11.multilocus import MultiLocusGeneticValue
from fwdpy11.trait_values import SlocusAdditiveTrait
from collections import namedtuple


def make_parser():
    parser = argparse.ArgumentParser(
        description="Process multi-locus outputs")

    parser.add_argument('--tarfile', '-t', type=str,
                        help=".tar file name containing replicates")
    parser.add_argument('--nprocs', '-p', type=int, default=1,
                        help="Number of processes to use")
    return parser

def get_params(infile):
    x = infile.split('.')
    m = [i for i,v in enumerate(x) if v.find('mu') != -1]
    x[m[0]]=x[m[0]].replace('mu','')
    mu = float(x[m[0]]+'.'+x[m[0]+1])
    o = [i for i,v in enumerate(x) if v.find('opt') != -1]
    x[o[0]]=x[o[0]].replace('opt','')
    opt = float(x[o[0]]+'.'+x[o[0]+1])
    return mu, opt

DATUM=namedtuple('DATUM',['repid','mu','opt','generation','segsmall','seglarge','flarge','fsmall'])

def process_replicate(argtuple):
    args, infile = argtuple
    mu, opt = get_params(infile)
    ghat=2.*np.sqrt(2.)*np.sqrt(mu)
    tf = tarfile.open(args.tarfile, 'r')
    ti = tf.getmember(infile)
    lzma_file = tf.extract(ti)
    stats = []
    with lzma.open(infile, 'rb') as f:
        while True:
            try:
                repid, pop = pickle.load(f)
                assert len(pop.mcounts)==len(pop.mutations), "fatal error"
                ss=0
                sl=0
                fs=0
                fl=0
                for i,j in zip(pop.mcounts, pop.mutations):
                    if i > 0 and j.neutral is False:
                        if i < 2*pop.N:
                            if np.abs(j.s) < ghat:
                                ss += 1
                            else:
                                sl += 1
                        else:
                            if np.abs(j.s) < ghat:
                                fs += 1
                            else:
                                fl += 1
                stats.append(DATUM(repid,mu,opt,pop.generation,ss,sl,fl,fs))
                    
                # Get total number of mutations of large vs small effect
                # that are segregating vs fixed.
                # Also record fixations
            except EOFError:
                os.remove(infile)
                return stats
            except:
                os.remove(infile)
                raise


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    if os.path.exists(args.tarfile) is False:
        raise RuntimeError(args.tarfile + " could not be found")

    dbname = args.tarfile.replace('.tar', '.numbers_of_variants.sqlite3')
    if os.path.exists(dbname):
        os.remove(dbname)
    tf = tarfile.open(args.tarfile, 'r')
    files = tf.getnames()
    tf.close()
    raw_args = [(args, i) for i in files] 
    # x=process_replicate(raw_args[0])
    # print(x)
    # sys.exit(0)
    with concurrent.futures.ProcessPoolExecutor(max_workers=32) as executor:
        futures = {executor.submit(process_replicate, i): i for i in raw_args}
        for fut in concurrent.futures.as_completed(futures):
            stats = fut.result()
            conn = sqlite3.connect(dbname)
            df = pd.DataFrame(stats)
            df.to_sql('data', conn, if_exists="append", index=False)
            conn.close()
