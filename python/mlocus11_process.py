import sys
import os
import pickle
import tarfile
import argparse
import lzma
import concurrent.futures
import numpy as np
import pandas as pd
from libsequence.polytable import SimData
from libsequence.windows import Windows
from libsequence.summstats import PolySIM,garudStats,std_nSLiHS
from libsequence.parallel import Scheduler
from fwdpy11.sampling import sample_separate

def make_parser():
    parser = argparse.ArgumentParser(
        description="Process multi-locus outputs")

    parser.add_argument('--tarfile','-t',type=str,help=".tar file name containing replicates")
    parser.add_argument('--nprocs','-p',type=int,default=1,help="Number of processes to use")
    parser.add_argument('--outfile','-o',type=str,help="Output file name.")
    parser.add_argument('--nsam','-n',type=int,default=50,help="Sample size (no. diploids)")
    parser.add_argument('--seed','-s',type=int,default=42,help="Random number seed")
    return parser

def make_empty_summstats_array():
    rv=np.array([(-1,-1,-1,-1,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)]*11,
            dtype=[('locus',np.int32),
                ('window',np.int32),
                ('repid',np.int32),
                ('generation',np.int32),
                ('tajd',np.float),
                ('pi',np.float),
                ('H1',np.float),
                ('H12',np.float),
                ('H2H1',np.float),
                ('nSL',np.float),
                ('iHS',np.float)])
    return rv

def get_summstats(pop,repid,nsam):
    ind = np.random.choice(pop.N,nsam,replace=False)
    s = sample_separate(pop,ind)
    temp=make_empty_summstats_array()
    rv=np.array([],temp.dtype)
    locus=0
    for si,bi in zip(s,pop.locus_boundaries[1:]):
        sd=SimData(si[0])
        w=Windows(sd,1.0,1.0,bi[0],bi[1])
        for i in range(len(w)):
            ps=PolySIM(w[i])
            gs=garudStats(w[i])
            nSL=std_nSLiHS(w[i],0.05,0.1)
            temp[i]=(locus,int(i),repid,pop.generation,ps.tajimasd(),ps.thetapi(),
                gs['H1'],gs['H12'],gs['H2H1'],nSL[0],nSL[1])
        rv=np.append(rv,temp.copy())
        locus+=1
    return rv


def process_replicate(argtuple):
    seed,args,infile=argtuple
    np.random.seed(seed)
    sched=Scheduler(1)
    tf=tarfile.open(args.tarfile,'r')
    ti=tf.getmember(infile)
    lzma_file=tf.extract(ti)
    summstats=[]
    with lzma.open(infile,'rb') as f:
        while True:
            try:
                rec = pickle.load(f)
                ss=get_summstats(rec[1],rec[0],args.nsam)
                summstats.append(ss)
                print(len(summstats))
                rec[1].clear()
            except:
                break
    
    return True

if __name__ == "__main__":
    parser=make_parser()
    args=parser.parse_args(sys.argv[1:])
    initial_seed = args.seed

    if os.path.exists(args.tarfile) is False:
        raise RuntimeError(args.tarfile + " could not be found")

    tf=tarfile.open(args.tarfile,'r')
    files = tf.getnames()
    tf.close()
    raw_args=[(i,args,j) for i,j in zip(np.random.choice(int(4e6),len(files),replace=False),files)]
    print(raw_args[0])
    rv=process_replicate(raw_args[0])
    print(rv)
    sys.exit(0)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.nprocs) as executor:
        futures = {executor.submit(process_replicate,i): i for i in raw_args[:1]}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            print(fn)
