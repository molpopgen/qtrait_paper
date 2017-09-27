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

def make_parser():
    parser = argparse.ArgumentParser(
        description="Process multi-locus outputs")

    parser.add_argument('--tarfile','-t',type=str,help=".tar file name containing replicates")
    parser.add_argument('--nprocs','-p',type=int,default=1,help="Number of processes to use")
    parser.add_argument('--outfile','-o',type=str,help="Output file name.")
    parser.add_argument('--nsam','-n',type=int,default=50,help="Sample size (no. diploids)")
    parser.add_argument('--seed','-s',type=int,default=42,help="Random number seed")
    return parser


def process_replicate(argtuple):
    args,infile=argtuple
    tf=tarfile.open(args.tarfile,'r')
    ti=tf.getmember(infile)
    lzma_file=tf.extract(ti)
    fixation = None
    with lzma.open(infile,'rb') as f:
        while True:
            try:
                repid, pop = pickle.load(f)
                fixations = np.array(pop.fixations.array())
            except:
                break
    os.remove(infile) 
    return (repid,fixations)

if __name__ == "__main__":
    parser=make_parser()
    args=parser.parse_args(sys.argv[1:])
    if os.path.exists(args.tarfile) is False:
        raise RuntimeError(args.tarfile + " could not be found")

    dbname = args.tarfile.replace('.tar','.fixations.db')
    #for d in dbnames:
    if os.path.exists(dbname):
        os.remove(dbname)
    tf=tarfile.open(args.tarfile,'r')
    files = tf.getnames()
    tf.close()
    raw_args=[(args,j) for j in files]
    #x=process_replicate(raw_args[0])
    #print(x)
    #sys.exit(0)
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        futures = {executor.submit(process_replicate,i): i for i in raw_args}
        for fut in concurrent.futures.as_completed(futures):
            repid, fixations = fut.result()
            conn = sqlite3.connect(dbname)
            df=pd.DataFrame(data)
            df.to_sql('data',conn,if_exists="append",index=False)
            conn.close()
