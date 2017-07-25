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
from libsequence.polytable import SimData
from libsequence.windows import Windows
from libsequence.summstats import PolySIM,garudStats,std_nSLiHS
from libsequence.parallel import Scheduler
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

def make_empty_summstats_array():
    rv=np.array([tuple([-1]*4+[np.nan]*8)]*11,
                dtype=[('locus',np.int32),
                ('window',np.int32),
                ('repid',np.int32),
                ('generation',np.int32),
                ('tajd',np.float),
                ('hprime',np.float),
                ('pi',np.float),
                ('H1',np.float),
                ('H12',np.float),
                ('H2H1',np.float),
                ('nSL',np.float),
                ('iHS',np.float)])
    return rv

def get_summstats(pop,repid,nsam,temp):
    """
    The 'genome scan' bit.
    """
    ind = np.random.choice(pop.N,nsam,replace=False)
    s = sample_separate(pop,ind)
    rv=np.array([],dtype=temp.dtype)
    locus=0
    for si,bi in zip(s,pop.locus_boundaries):
        sd=SimData(si[0])
        w=Windows(sd,1.0,1.0,bi[0],bi[1])
        for i in range(len(w)):
            ps=PolySIM(w[i])
            gs=garudStats(w[i])
            nSL=std_nSLiHS(w[i],0.05,0.1)
            temp[i]=np.array([(locus,int(i),repid,pop.generation,ps.tajimasd(),
                    ps.thetapi(),ps.hprime(),
                    gs['H1'],gs['H12'],gs['H2H1'],nSL[0],nSL[1])],dtype=temp.dtype)
            #print("done assigning ",i)
        #print("attempt concat")
        rv=np.concatenate([rv,temp.copy()])
        #print("here at ",locus)
        locus+=1
    #print("returning summstats")
    return rv


def get_genetic_values_per_locus(gv_fxn,pop,temp):
    """
    Return the mean and variance of G, separated
    per locus and the VG

    This is not equivalent to the VG due to a locus,
    b/c E[sum(x)] != sum(E[x]).
    """
    j=0
    for i in pop.diploids:
        temp[j] = gv_fxn(i,pop)
        j+=1
    csums = temp.sum(1)
    return(temp.mean(0),csums.var())

def process_replicate(argtuple):
    seed,args,infile=argtuple
    np.random.seed(seed)
    sched=Scheduler(1)
    tf=tarfile.open(args.tarfile,'r')
    ti=tf.getmember(infile)
    lzma_file=tf.extract(ti)
    gvtemp = None
    vg_over_time = np.array([],dtype=[('repid',np.int32),
        ('generation',np.int32),
        ('vg',np.float)])
    gvalues = np.array([],dtype=[('repid',np.int32),('generation',np.int32),('locus',np.int32),('g',np.float)])
    gvalues_intermediate = np.array([(-1,-1,-1,np.nan)]*10,dtype=gvalues.dtype)
    genetic_value_fxn = MultiLocusGeneticValue([SlocusAdditiveTrait(2.0)]*10)
    genome_scan_temp=make_empty_summstats_array()
    summstats=np.array([],dtype=genome_scan_temp.dtype)
    with lzma.open(infile,'rb') as f:
        while True:
            try:
                rec = pickle.load(f)
                if gvtemp is None or gvtemp.shape[0] != rec[1].N:
                    gvtemp = np.zeros((rec[1].N,len(rec[1].locus_boundaries)))
                #print("here")
                gv_vg = get_genetic_values_per_locus(genetic_value_fxn,
                        rec[1],gvtemp)
                #print("here2")
                #print(gv_vg)
                #print(gv_vg[1])
                #print(rec[0],rec[1].generation,gv_vg[1])
                cc=np.array([(rec[0],rec[1].generation,gv_vg[1])])
                #print(cc.dtype)
                #print(vg_over_time.dtype)
                vg_over_time = np.concatenate([vg_over_time,
                    np.array([(rec[0],rec[1].generation,gv_vg[1])],dtype=vg_over_time.dtype)])
                #print("here2b")
                for i in range(len(gv_vg[0])):
                    gvalues_intermediate[i]=(rec[0],rec[1].generation,i,gv_vg[0][i])
                #print("here3")
                gvalues=np.concatenate([gvalues,gvalues_intermediate])
                #print("here4")
                ss=get_summstats(rec[1],rec[0],args.nsam,genome_scan_temp)
                #print(ss)
                #print(ss.dtype)
                #print(summstats.dtype)
                summstats=np.concatenate([summstats,ss])
                #print("here5")
                rec[1].clear()
            except:
                break
    os.remove(infile) 
    return (summstats,gvalues,vg_over_time)

if __name__ == "__main__":
    parser=make_parser()
    args=parser.parse_args(sys.argv[1:])
    np.random.seed(args.seed)
    if os.path.exists(args.tarfile) is False:
        raise RuntimeError(args.tarfile + " could not be found")

    dbnames = [args.tarfile.replace('.tar','.genome_scan.db'),
            args.tarfile.replace('.tar','.genetic_values_per_locus.db'),
            args.tarfile.replace('.tar','.vg_over_time.db')]
    #for d in dbnames:
    #    if os.path.exists(d):
    #        os.remove(d)
    tf=tarfile.open(args.tarfile,'r')
    files = tf.getnames()
    tf.close()
    raw_args=[(i,args,j) for i,j in zip(sorted(np.random.choice(int(4e6),len(files),replace=False)),files)]
    #x=process_replicate(raw_args[0])
    #print(x)
    #sys.exit(0)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.nprocs) as executor:
        futures = {executor.submit(process_replicate,i): i for i in raw_args}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            for dbname,data in zip(dbnames,fn):
                conn = sqlite3.connect(dbname)
                df=pd.DataFrame(data)
                df.to_sql('data',conn,if_exists="append",index=False)
                conn.close()
    for dbname in dbnames:
        conn = sqlite3.connect(dbname)
        conn.execute("create index gen_rep if not exists on data")
        conn.close()
