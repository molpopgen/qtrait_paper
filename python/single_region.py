from __future__ import print_function
import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import gzip,argparse,sys,warnings,math
import pandas as pd
##import local functions
from single_region_common import *
import pyximport
pyximport.install()
import loads_opt
SAMPLERS=('stats','freq','load')
TRAITS=('additive','mult')
trait_models = {'additive':qt.SpopAdditiveTrait(),
                'mult':qt.SpopMultTrait()}

def get_trait_model(traitString):
   return trait_models[traitString]

def get_sampler_type(samplerString, traitString,length,optimum):
    if samplerString == 'stats':
        return fp.QtraitStatsSampler(length,optimum)
    elif samplerString == 'freq':
        return fp.FreqSampler(length)
    elif samplerString == 'load':
        if traitString == 'additive':
            return loads_opt.additiveLoad(length,optimum)
        elif traitString == 'mult':
            return loads_opt.multiplicativeLoad(length,optimum)
    else:
        print ("invalid sampler name")
        usage()
        sys.exit(0)

def write_output(sampler,outputFilename,REPID,batch,mode,dfname):
    if isinstance(sampler,fp.FreqSampler):
        ##Create a new file for each replicate so that file sizes
        ##don't get unwieldy
        fn=outputFilename+'.batch'+str(batch)+'.h5'
	output = pd.HDFStore(fn,mode,complevel=6,complib='zlib')	
        for di in sampler:
            df=pd.DataFrame(fp.tidy_trajectories(di))
            df['rep']=REPID*len(df.index)
            df.reset_index(['rep'],inplace=True,drop=True)
            output.append(dfname,df)
            REPID+=1
        output.close()
    else:
        #elif isinstance(sampler,fp.QtraitStatsSampler):
        ##Write in append mode
        output = pd.HDFStore(outputFilename,mode,complevel=6,complib='zlib')
        for s in sampler:
            df=pd.DataFrame(s)
            df['rep']=[REPID]*len(df.index)
            REPID+=1
            output.append(dfname,df)
        output.close()
    #else:
    #    raise RuntimeError("uh oh: sampler type not recognized for output.  We shouldn't have gotten this far!")
    return REPID
       
def make_parser():
    parser=argparse.ArgumentParser(description="Single region, simple demography, additive effects.")

    parser.add_argument('--verbose',action='store_true')
    group = parser.add_argument_group("Required arguments")
    group.add_argument('--mutrate','-m',type=float,default=-1.0,help="Mutation rate to selected variants. Required.",required=True)
    group.add_argument("--outfile",'-O',type=str,default=None,help="Output file. Required",required=True)
    group.add_argument('--sampler',choices=SAMPLERS,help="Temporal sampler type.",required=True)

    parser.add_argument('--seed','-S',type=int,default=0,metavar='SEED',help='Random number seed (default 0)')
    #parser.add_argument('--theta','-th',type=float,default=100.,metavar='THETA',help="4Nu/region (default 100.0)")
    parser.add_argument('--recrate','-r',type=float,default=0.5,help="Genetic distance of region.  Default = 0.5 = 50cM")
    parser.add_argument('--tsample','-t',type=int,default=0,metavar='TSAMPLE',help="How often to apply sampler. Default = 0 = never")
    parser.add_argument('--sigmu',type=float,default=0.25,metavar='SIGMU',help="Std. dev. of effect sizes. Default = 0.25")
    parser.add_argument('--H2',type=float,default=1.0,help="Desired heritability.  Default of 1.0. We use HoC approximation to determine std. dev. of environmental/random effects to give this heritability on average")
    parser.add_argument('--VS',type=float,default=1.0,help="Total variation in selection intensity (strenth of selection against variance in trait value), or V(S). Default to 1.0. If H2=1, then V(S) is interpreted as V(S)+V(E).")
    parser.add_argument('--trait',choices=TRAITS,default='additive',help="Genotype to phenotype model. Default is additive.")
    parser.add_argument('--dominance',type=float,default=1.0,help="Dominance. Default is 1.0 (co-dominance).")
    parser.add_argument('--ncores',type=int,default=64,help="No. simultaneous simulations to run. Default = 64")
    parser.add_argument('--nbatches',type=int,default=16,help="No. 'batches' of simulations to run.  Total number of replicates will be NCORES*NBATCHES")
    #Positional arguments
    parser.add_argument('popsize',type=int,help="Population size.  Simulation will run 10*N generations before and after optimum shift.")
    parser.add_argument('optimum',type=float,help="Optimum trait value.")

    return parser

def main():
    parser=make_parser()
    args=parser.parse_args(sys.argv[1:])

    if args.verbose:
        print (args)
    ##Figure out sigma_E from params
    sigE = get_sigE_additive(args.mutrate,args.VS,args.H2)

    trait = get_trait_model(args.trait)

    nlist = np.array([args.popsize]*(10*args.popsize),dtype=np.uint32)
    rng=fp.GSLrng(args.seed)
    mu_neutral = 0.0
    nregions=[]
    sregions=[fp.GaussianS(0,1,1,args.sigmu,args.dominance)]
    recregions=[fp.Region(0,1,1)]
    
    if args.sampler == 'stats':
        output=pd.HDFStore(args.outfile,'w',complevel=6,complib='zlib')
        output.close()
    REPID=0
    for BATCH in range(args.nbatches):
        pops=fp.SpopVec(args.ncores,args.popsize)
        sampler=get_sampler_type(args.sampler,args.trait,len(pops),0.0)
        qt.evolve_regions_qtrait_sampler_fitness(rng,pops,sampler,trait,
                                                 nlist,
                                                 0.0,
                                                 args.mutrate,
                                                 args.recrate,
                                                 nregions,
                                                 sregions,
                                                 recregions,
                                                 args.tsample,
                                                 sigE,
                                                 optimum=0.0,
                                                 VS=args.VS)
        if args.sampler == 'freq':
            dummy=write_output(sampler,args.outfile,REPID,BATCH,'w','trajectories')
        elif args.sampler == 'stats':
            dummy=write_output(sampler,args.outfile,REPID,BATCH,'a','stats')
        else:
            dummy=write_output(sampler,args.outfile,REPID,BATCH,'a','load')
        sampler=get_sampler_type(args.sampler,args.trait,len(pops),args.optimum)
        qt.evolve_regions_qtrait_sampler_fitness(rng,pops,sampler,trait,
                                                 nlist,
                                                 0.0,
                                                 args.mutrate,
                                                 args.recrate,
                                                 nregions,
                                                 sregions,
                                                 recregions,
                                                 args.tsample,
                                                 sigE,
                                                 optimum=args.optimum,
                                                 VS=args.VS)
        if args.sampler == 'freq':
            #Append this time!
            REPID=write_output(sampler,args.outfile,REPID,BATCH,'a','trajectories')
        elif args.sampler == 'stats':
            REPID=write_output(sampler,args.outfile,REPID,BATCH,'a','stats')
        else:
            REPID=write_output(sampler,args.outfile,REPID,BATCH,'a','load')

if __name__ == "__main__":
    main()
