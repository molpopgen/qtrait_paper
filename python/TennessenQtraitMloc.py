#Simulate the Tennessen model with the following params:
#1. 10 loci separated by 50cM
#2. Each locus is 11*theta and 11*rho long, where theta & rho are per-locus input params
#3. The middle of each locus is where causative mutations happen at rate mu with 
#   and sigmamu=0.25 and VS=1
#4. A sample of n=6e3 is taken at the end and summary stats are calcuated for 11 window
#5. The samples are written to files for latter analysis with SweeD and SDS
#6. The change in the optimum trait value occurs "out of Africa".  The optimum is 0 prior to then.
from __future__ import print_function
import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import fwdpy.demography as fpd
import numpy as np
import pandas as pd
import sys,gc
import datetime
import warnings
import argparse
import random
import multiprocessing as mp
import libsequence.parallel as lsp
import libsequence.polytable as pt
import libsequence.summstats as sstats
import libsequence.windows as windows
#do check on minimum pylibseq version
import libsequence
assert libsequence.__version__ >= '0.1.9'

STATNAMES=['S','tajimasd','pi','thetaw','hprime','H1','H12','H2H1','nSL','iHS']

def get_nlist1():
    """
    Generates a numpy array of the canges in N over time

    There are 5 epochs, with t=0 being the present.
    
    E1: Ne= 7,310  from t=start(8N?) to t = - 5920 generation (Ancestral sizes until 5920 generations ago)
    
    E2: Ne =14,474 from t = -5920 to t = -2040 (Ancient growth at 5920 g ago)
    
    E3: Ne =1,861 from t= -2040 to t= -920 (OOA, first bottle neck 2040 g ago)
    
    E4: Ne = 1,032 to Ne = 9,300 during t = -920 to t = -205 ( second bottle neck and onset of 715 g of exponential growth at rate 0.31% per gen )  
    
    E5: Ne = 9,300 to Ne = 512,000 during t = -205 to t = -0 ( 205 g of exponential growth at rate 1.95% per gen )  
    """
    n=[7310]*(10*7310) #E1: evolve ancestral size to mutation/selection/drift equilibrium
    n.extend([14474]*(5920-2040)) #E2
    #n.extend([1861]*(2040-920)) #E3
    #n.extend(fpd.exponential_size_change(1032,9300,920-205)) #E4
    #n.extend(fpd.exponential_size_change(9300,51200,205)) #E5
    return n

def get_nlist2():
    """
    Generates a numpy array of the canges in N over time

    There are 5 epochs, with t=0 being the present.
    
    E1: Ne= 7,310  from t=start(8N?) to t = - 5920 generation (Ancestral sizes until 5920 generations ago)
    
    E2: Ne =14,474 from t = -5920 to t = -2040 (Ancient growth at 5920 g ago)
    
    E3: Ne =1,861 from t= -2040 to t= -920 (OOA, first bottle neck 2040 g ago)
    
    E4: Ne = 1,032 to Ne = 9,300 during t = -920 to t = -205 ( second bottle neck and onset of 715 g of exponential growth at rate 0.31% per gen )  
    
    E5: Ne = 9,300 to Ne = 512,000 during t = -205 to t = -0 ( 205 g of exponential growth at rate 1.95% per gen )  
    """
    n=[1861]*(2040-920) #E3
    n.extend(fpd.exponential_size_change(1032,9300,920-205)) #E4
    n.extend(fpd.exponential_size_change(9300,512000,205)) #E5
    return n

def process_samples(args):
    di,statfile,locus_boundaries,repid=args
    H5out = pd.HDFStore(statfile,mode='a'i,complevel=6,complib='zlib')
    #get summary stats in sliding windows based on neutral diversity
    sd=[pt.SimData(dii[0][0]) for dii in di]
    w=[windows.Windows(i,window_size=1.,step_len=1.,starting_pos=j[0],ending_pos=j[1]) for i,j in zip(sd,locus_boundaries)]
    sd=None
    hapstats=[] #nSL, etc.
    polySIMlist=[]
    for wi in w:
        polySIMlist.extend([sstats.PolySIM(j) for j in wi])
        hapstats.extend([(sstats.garudStats(j),sstats.std_nSLiHS(j,0.05,0.1)) for j in wi])
    wi = None #Delete the last window...
    #stats=[(i.numpoly(),i.tajimasd(),i.thetapi(),i.thetaw(),i.hprime()) for i in polySIMlist]
    stats=[(i.numpoly(),i.tajimasd(),i.thetapi(),i.thetaw(),i.hprime()) for i in polySIMlist]
    combinedStats=[i+(j[0]['H1'],j[0]['H12'],j[0]['H2H1'],j[1][0],j[1][1]) for i,j in zip(stats,hapstats)]
    polySIMlist=None
    #combinedStats=[i+(j[0]['H1'],j[0]['H12'],j[0]['H2H1'],j[1][0],j[1][1]) for i,j in zip(stats,hapstats)]
    window=0
    locus=0
    for cs in combinedStats:
        tempDF=pd.DataFrame({'stat':STATNAMES,'value':list(cs),
                'window':[window]*len(STATNAMES),'rep':[repid]*len(STATNAMES),
                'locus':[locus]*len(STATNAMES)})
                #'locus':[l for l in range(10) for i in range(len(locus_boundaries))]})
        H5out.append('stats',tempDF)
        window+=1
        if window==11:
            window=0
            locus+=1
        tempDF=None
    cs=None
    combinedStats=None
    stats=None
    hapstats=None
    H5out.close()

def make_parser():
    parser = argparse.ArgumentParser(
        description="Ten locus simulation of Tennessen model")

    parser.add_argument('--mu','-m',type=float,default=1e-3,help="mutation rate to mutations affecting trait value")
    parser.add_argument('--theta','-t',type=float,default=100.0,metavar="THETA",help="4Nu")
    parser.add_argument('--rho','-r',type=float,default=100.0,metavar="RHO",help="4Nr")
    parser.add_argument('--TBB','-T',type=int,default=-1,help="TBB threading control")
    parser.add_argument('--seed','-S',type=int,default=None,help="RNG seed")
    parser.add_argument('--ncores','-n',type=int,default=1)
    parser.add_argument('--nreps','-nr',type=int,default=1)
    parser.add_argument('--nloci','-nl',type=int,default=10)
    parser.add_argument('--sigmu','-si',type=float,default=0.25)
    parser.add_argument('--statfile','-st',type=str,default=None)
    parser.add_argument('--opt','-o',type=float,default=0.0)

    return parser

def run_batch(argtuple):
    args,repid,batch=argtuple
    print ("seed for batch = ",args.seed)
    nstub = "neutral.mu"+str(args.mu)+".opt"+str(args.opt)
    sstub = "selected.mu"+str(args.mu)+".opt"+str(args.opt)
    rnge=fp.GSLrng(args.seed)

    NANC=7310
    locus_boundaries=[(float(i+i*11),float(i+i*11+11)) for i in range(args.nloci)]
    nregions=[fp.Region(j[0],j[1],args.theta/(4.*float(NANC)),coupled=True) for i,j in zip(range(args.nloci),locus_boundaries)]
    recregions=[fp.Region(j[0],j[1],args.rho/(4.*float(NANC)),coupled=True) for i,j in zip(range(args.nloci),locus_boundaries)]
    sregions=[fp.GaussianS(j[0]+5.,j[0]+6.,args.mu,args.sigmu,coupled=False) for i,j in zip(range(args.nloci),locus_boundaries)]
    f=qtm.MlocusAdditiveTrait()

    nlist=np.array(get_nlist1(),dtype=np.uint32)
    pops = fp.MlocusPopVec(args.ncores,nlist[0],args.nloci)
    sampler=fp.NothingSampler(len(pops))
    d=datetime.datetime.now()
    print("starting batch, ",batch, "at ",d.now())
    qtm.evolve_qtraits_mloc_regions_sample_fitness(rnge,pops,sampler,f,
            nlist[0:],
            nregions,sregions,recregions,
            [0.5]*(args.nloci-1),
            0,
            0,0.)
    d=datetime.datetime.now()
    print(d.now())
    nlist=np.array(get_nlist2(),dtype=np.uint32)
    qtm.evolve_qtraits_mloc_regions_sample_fitness(rnge,pops,sampler,f,
            nlist[0:],
            nregions,sregions,recregions,
            [0.5]*(args.nloci-1),
            0,args.opt)
    d=datetime.datetime.now()
    print(d.now())
    if args.statfile is not None:
        sched = lsp.scheduler_init(args.TBB)
    for pi in pops:
        neutralFile = nstub + '.rep' + str(repid) + '.gz'
        selectedFile = sstub + '.rep' + str(repid) + '.gz'
        BIGsampler=fp.PopSampler(1,6000,rnge,False,neutralFile,selectedFile,recordSamples=True,boundaries=locus_boundaries)
        fp.apply_sampler_single(pi,BIGsampler)
        if args.statfile is not None:
            for di in BIGsampler:
                process_samples((di,args.statfile,locus_boundaries,repid))
        repid+=1
    #neutralFile = nstub + '.batch' + str(batch)
    #selectedFile = sstub + '.batch' + str(batch)
    #BIGsampler=fp.PopSampler(len(pops),600,rnge,False,neutralFile,selectedFile,recordSamples=True,boundaries=locus_boundaries)
    #fp.apply_sampler(pops,BIGsampler)
    #d=datetime.datetime.now()
    #print(d.now())
    #if args.statfile is not None:
    #    sched = lsp.scheduler_init(args.TBB)
    #    for di in BIGsampler:
    #        process_samples((di,args.statfile,locus_boundaries,repid))
    #        repid+=1
    #        del di
    #BIGsampler.force_clear()
    #BIGsampler=None
    pops.clear()
    pops=None

if __name__ == "__main__":
    parser=make_parser()
    args=parser.parse_args(sys.argv[1:])
    if args.statfile is not None:
        H5out = pd.HDFStore(args.statfile,'w',complevel=6,complib='zlib')
        H5out.close()
    #Use multiprocessing hack to keep total RAM well-controlled.
    repid=0
    #Use user seed to generate a 
    #seed for each "batch" of "ncore" reps
    initial_seed = args.seed
    random.seed(initial_seed)
    for batch in range(args.nreps):
        repseed=random.randrange(42000000)
        args.seed=repseed
        P=mp.Pool(1)
        P.map(run_batch,[(args,repid,batch)])
        P.close()
        repid+=args.ncores
