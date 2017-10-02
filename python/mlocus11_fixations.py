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
    return parser

def make_df(data):
    repid, fixations, fixation_times = data
    dft = pd.DataFrame(fixations)
    dft.loc[:,'ftime'] = fixation_times
    df = dft.query('neutral == 0')
    arose_after = df.g > 50000
    arose_before = df.g < 50000
    fixed_after = df.ftime > 50000
    df.loc[:,'sweep_type'] = 'none'
    df.loc[arose_after,'sweep_type'] = 'hard'
    df.loc[arose_before&fixed_after,'sweep_type'] = 'soft'
    df.loc[:,'repid'] = [repid]*len(df.index)
    df.reset_index(inplace=True)
    df.drop(['index','label','neutral'],inplace=True,axis=1)
    return df

def process_replicate(argtuple):
    args,infile=argtuple
    tf=tarfile.open(args.tarfile,'r')
    ti=tf.getmember(infile)
    lzma_file=tf.extract(ti)
    fixations = None
    with lzma.open(infile,'rb') as f:
        while True:
            try:
                repid, pop = pickle.load(f)
                fixations = np.array(pop.fixations.array())
                fixation_times = pop.fixation_times
                # if any(i.neutral == 0 for i in pop.fixations) is True:
                #     print("doit!")
                #     print(make_df((repid,fixations,fixation_times)))
                #     print("done")
                # zbar = np.array(pop.diploids.trait_array())['g'].mean()
                # print(pop.generation,zbar)
            except:
                break
    os.remove(infile) 
    return make_df((repid, fixations, fixation_times))

if __name__ == "__main__":
    parser=make_parser()
    args=parser.parse_args(sys.argv[1:])
    if os.path.exists(args.tarfile) is False:
        raise RuntimeError(args.tarfile + " could not be found")

    dbname = args.tarfile.replace('.tar','.fixations.db')
    if os.path.exists(dbname):
        os.remove(dbname)
    tf=tarfile.open(args.tarfile,'r')
    files = tf.getnames()
    tf.close()
    raw_args=[(args,j) for j in files]
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        futures = {executor.submit(process_replicate,i): i for i in raw_args}
        for fut in concurrent.futures.as_completed(futures):
            df = fut.result()
            conn = sqlite3.connect(dbname)
            df.to_sql('data',conn,if_exists="append",index=False)
            conn.close()
