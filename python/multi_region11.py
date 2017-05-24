#1. 10 loci separated by 50cM
#2. Each locus is 11*theta and 11*rho long, where theta & rho are per-locus input params
#3. The middle of each locus is where causative mutations happen at rate mu with 
#   and sigmamu=0.25 and VS=1
import fwdpy11 as fp11
import fwdpy11.wright_fisher_qtrait as fp11qt
import fwdpy11.model_params as fp11mp
import fwdpy11.trait_values as fp11tv
import fwdpy11.multilocus as fp11ml
import numpy as np
import pandas as pd
import sys,gc
import datetime
import warnings
import argparse
import random
import pickle
import os
import lzma
import tarfile
import concurrent.futures
#import multiprocessing as mp

def make_parser():
    parser = argparse.ArgumentParser(
        description="Simple multi-locus simulation of single optimum shift")

    parser.add_argument('--N','-N',type=int,default=5000,help="Diploid population size")
    parser.add_argument('--mu','-m',type=float,default=1e-3,help="mutation rate to mutations affecting trait value")
    parser.add_argument('--theta','-t',type=float,default=100.0,metavar="THETA",help="4Nu")
    parser.add_argument('--rho','-r',type=float,default=100.0,metavar="RHO",help="4Nr")
    parser.add_argument('--seed','-S',type=int,default=None,help="RNG seed")
    parser.add_argument('--ncores','-n',type=int,default=1)
    parser.add_argument('--nreps','-nr',type=int,default=1)
    parser.add_argument('--nloci','-nl',type=int,default=10)
    parser.add_argument('--sigmu','-si',type=float,default=0.25)
    parser.add_argument('--opt','-o',type=float,default=0.0)
    parser.add_argument('--stub','-s',type=str)
    parser.add_argument('--tarfile',type=str,default=None,help="If not None, merge all output into a single tar file.")
    return parser

class Pickler(object):
    def __init__(self,repid,N,filename):
        self.repid=repid
        self.N=N
        self.filename=filename
        self.last_gen_recorded=-1
        if os.path.exists(self.filename):
            os.remove(self.filename)
    def __call__(self,pop):
        c1=(pop.generation >= 8*self.N and pop.generation < 10*self.N and pop.generation % 500 == 0.0)
        c2=(pop.generation >= 10*self.N and pop.generation <= 12*self.N and pop.generation % 10 == 0.0)
        c3=(pop.generation > 12*self.N and pop.generation <= 14*self.N and pop.generation % 50 == 0.0)
        c4=(pop.generation > 14*self.N and pop.generation % 500 == 0.0)
        if c1 or c2 or c3 or c4: 
            with lzma.open(self.filename,"ab",preset=9) as f:
                pickle.dump((self.repid,pop),f,-1)
            self.last_gen_recorded = pop.generation

def run_replicate(argtuple):
    args,repid,repseed=argtuple
    rng=fp11.GSLrng(repseed)
    NANC=args.N
    locus_boundaries=[(float(i+i*11),float(i+i*11+11)) for i in range(args.nloci)]
    nregions=[[fp11.Region(j[0],j[1],args.theta/(4.*float(NANC)),coupled=True)] for i,j in zip(range(args.nloci),locus_boundaries)]
    recregions=[[fp11.Region(j[0],j[1],args.rho/(4.*float(NANC)),coupled=True)] for i,j in zip(range(args.nloci),locus_boundaries)]
    sregions=[[fp11.GaussianS(j[0]+5.,j[0]+6.,args.mu,args.sigmu,coupled=False)] for i,j in zip(range(args.nloci),locus_boundaries)]
    interlocus_rec=fp11ml.binomial_rec(rng,[0.5]*(args.nloci-1))
    nlist=np.array([NANC]*20*NANC,dtype=np.uint32)
    env = [(0,0,1),(10*NANC,args.opt,1)]
    pdict = {'nregions':nregions,
            'sregions':sregions,
            'recregions':recregions,
            'demography':nlist,
            'interlocus':interlocus_rec,
            'agg':fp11ml.AggAddTrait(),
            'gvalue':fp11ml.MultiLocusGeneticValue([fp11tv.SlocusAdditiveTrait(2.0)]*args.nloci),
            'trait2w':fp11qt.GSSmo(env),
            'mutrates_s':[args.mu/float(args.nloci)]*args.nloci,
            'mutrates_n':[float(args.nloci)*args.theta/float(4*NANC)]*args.nloci,
            'recrates':[float(args.nloci)*args.rho/float(4*NANC)]*args.nloci,
            }
    params=fp11.model_params.MlocusParamsQ(**pdict)
    ofilename = args.stub + '.rep' + str(repid) + '.pickle.lzma'
    recorder=Pickler(repid,NANC,ofilename)
    pop=fp11.MlocusPop(NANC,args.nloci,locus_boundaries)
    fp11qt.evolve(rng,pop,params,recorder)
    #Make sure last gen got pickled!
    if recorder.last_gen_recorded != pop.generation:
        with open(ofilename,"ab") as f:
            pickle.dump((repid,pop),f,-1)
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
    #run_replicate(arglist[0])
    #P = mp.Pool(args.ncores)
    #res = P.imap_unordered(run_replicate,arglist)
    #P.close()
    #P.join()
