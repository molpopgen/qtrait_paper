#Simulate the Tennessen model with the following params:
#1. 10 loci separated by 50cM
#2. Each locus is 11*theta and 11*rho long, where theta & rho are per-locus input params
#3. The middle of each locus is where causative mutations happen at rate mu with 
#   and sigmamu=0.25 and VS=1
#4. The change in the optimum trait value occurs "out of Africa".  The optimum is 0 prior to then.
import fwdpy11 as fp11
import fwdpy11.wright_fisher_qtrait as fp11qt
import fwdpy11.demography as fp11d
import fwdpy11.model_params as fp11mp
import fwdpy11.trait_values as fp11tv
import fwdpy11.multilocus as fp11ml
import fwdpy11.sampling
import numpy as np
import pandas as pd
import sys,gc
import datetime
import warnings
import argparse
import random
import pickle
import os
import gzip
import tarfile
import concurrent.futures
#import multiprocessing as mp

def get_nlist():
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
    n.extend([1861]*(2040-920)) #E3
    n.extend(fp11d.exponential_size_change(1032,9300,920-205)) #E4
    n.extend(fp11d.exponential_size_change(9300,512000,205)) #E5
    return np.array(n,dtype=np.uint32)

def make_parser():
    parser = argparse.ArgumentParser(
        description="Ten locus simulation of Tennessen model")

    parser.add_argument('--mu','-m',type=float,default=1e-3,help="mutation rate to mutations affecting trait value")
    parser.add_argument('--theta','-t',type=float,default=100.0,metavar="THETA",help="4Nu")
    parser.add_argument('--rho','-r',type=float,default=100.0,metavar="RHO",help="4Nr")
    parser.add_argument('--seed','-S',type=int,default=None,help="RNG seed")
    parser.add_argument('--nsam',type=int,default=250,help="Number of diploids to sample.")
    parser.add_argument('--ncores','-n',type=int,default=1)
    parser.add_argument('--nreps','-nr',type=int,default=1)
    parser.add_argument('--nloci','-nl',type=int,default=10)
    parser.add_argument('--sigmu','-si',type=float,default=0.25)
    parser.add_argument('--opt','-o',type=float,default=0.0)
    parser.add_argument('--stub','-s',type=str)
    parser.add_argument('--tarfile',type=str,default=None,help="If not None, merge all output into a single tar file.")
    return parser

def get_matrix(pop,ndips):
    dips = np.random.choice(pop.N,ndips,replace=False)
    keys = fwdpy11.sampling.mutation_keys(pop,dips)
    neutral_sorted_keys=[i for i in sorted(keys[0],key=lambda x,m=pop.mutations: m[x[0]].pos)]
    selected_sorted_keys=[i for i in sorted(keys[1],key=lambda x,m=pop.mutations: m[x[0]].pos)]
    m = fwdpy11.sampling.genotype_matrix(pop,dips,neutral_sorted_keys,selected_sorted_keys)
    gw = [(pop.diploids[i][0].g,pop.diploids[i][0].w) for i in dips] 
    return (m,gw)

class Pickler(object):
    def __init__(self,repid,nsam,gen1,gen2,filename):
        self.nsam=nsam
        self.repid=repid
        self.gen1=gen1
        self.gen2=gen2
        self.filename=filename
        self.last_gen_recorded=-1
        if os.path.exists(self.filename):
            os.remove(self.filename)
    def __call__(self,pop):
        if (pop.generation >= self.gen1 and pop.generation < self.gen2 and pop.generation % 500 == 0.0) or (pop.generation >= self.gen2 and pop.generation % 50 == 0.0):
            with gzip.open(self.filename,"ab") as f:
                m=get_matrix(pop,self.nsam)
                pickle.dump(m,f,-1)
                #pickle.dump((self.repid,pop),f,-1)
            self.last_gen_recorded = pop.generation

def run_replicate(argtuple):
    args,repid,repseed=argtuple
    rng=fp11.GSLrng(repseed)
    NANC=7310
    locus_boundaries=[(float(i+i*11),float(i+i*11+11)) for i in range(args.nloci)]
    nregions=[[fp11.Region(j[0],j[1],args.theta/(4.*float(NANC)),coupled=True)] for i,j in zip(range(args.nloci),locus_boundaries)]
    recregions=[[fp11.Region(j[0],j[1],args.rho/(4.*float(NANC)),coupled=True)] for i,j in zip(range(args.nloci),locus_boundaries)]
    sregions=[[fp11.GaussianS(j[0]+5.,j[0]+6.,args.mu,args.sigmu,coupled=False)] for i,j in zip(range(args.nloci),locus_boundaries)]
    interlocus_rec=fp11ml.binomial_rec(rng,[0.5]*(args.nloci-1))
    nlist=get_nlist()
    env = [(0,0,1),(np.argmax(nlist==1861),args.opt,1)]
    pdict = {'nregions':nregions,
            'sregions':sregions,
            'recregions':recregions,
            'demography':nlist,
            'interlocus':interlocus_rec,
            'agg':fp11ml.AggAddTrait(),
            'gvalue':fp11ml.MultiLocusGeneticValue([fp11tv.SlocusAdditiveTrait(2.0)]*args.nloci),
            'trait2w':fp11qt.GSSmo(env),
            'mutrates_s':[args.mu]*args.nloci,
            'mutrates_n':[float(args.nloci)*args.theta/float(4*NANC)]*args.nloci,
            'recrates':[float(args.nloci)*args.rho/float(4*NANC)]*args.nloci,
            }
    #print(pdict['mutrates_n'])
    #for i in sregions:
    #    for j in i:
    #        print(str(j))
    #return
    params=fp11.model_params.MlocusParamsQ(**pdict)
    ofilename = args.stub + '.rep' + str(repid) + '.gz'
    recorder=Pickler(repid,args.nsam,8*NANC,np.argmax(nlist == 1861),ofilename)
    pop=fp11.MlocusPop(NANC,args.nloci,locus_boundaries)
    fp11qt.evolve(rng,pop,params,recorder)
    #Make sure last gen got pickled!
    if recorder.last_gen_recorded != pop.generation:
        with open(ofilename,"ab") as f:
            with gzip.open(ofilename,"ab") as f:
                m=get_matrix(pop,args.nsam)
                pickle.dump(m,f,-1)
    pop.clear()
    del pop
    pop = None
    return ofilename

if __name__ == "__main__":
    parser=make_parser()
    args=parser.parse_args(sys.argv[1:])
    initial_seed = args.seed
    np.random.seed(initial_seed)
    repseeds = np.random.randint(0,42000000,args.nreps)
    arglist=[(args,i,j) for i,j in zip(range(len(repseeds)),repseeds)]
    tfn=args.tarfile
    if tfn is not None:
        tf=tarfile.open(tfn,"w")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.ncores) as executor:
        futures = {executor.submit(run_replicate,i): i for i in arglist}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            if tfn is not None:
                tf.add(fn)
                os.remove(fn)
    if tfn is not None:
        tf.close()
    #P = mp.Pool(args.ncores)
    #res = P.imap_unordered(run_replicate,arglist)
    #P.close()
    #P.join()
