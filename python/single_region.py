import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import gzip,argparse,getopt,sys,warnings,math
import pandas as pd
##import local functions
from single_region_common import *

SAMPLERS=('stats','freq')
TRAITS=('additive','mult')
#valid_trait_models=['additive','mult']
trait_models = {'additive':qt.SpopAdditiveTrait(),
                'mult':qt.SpopMultTrait()}

#valid_sampler_names=['stats','freq']
def usage():
    print "something went wrong"

def get_trait_model(traitString):
   return trait_models[traitString]

def get_sampler_type(samplerString, length,optimum):
    if samplerString == 'stats':
        return fp.QtraitStatsSampler(length,optimum)
    elif samplerString == 'freq':
        return fp.FreqSampler(length)
    else:
        print ("invalid sampler name")
        usage()
        sys.exit(0)

def write_output(sampler,outputFilename,REPID,batch,mode):
    if isinstance(sampler,fp.FreqSampler):
        ##Create a new file for each replicate so that file sizes
        ##don't get unwieldy
        fn=outputFilename+'.batch'+str(batch)+'.h5'
	output = pd.HDFStore(fn,mode,complevel=6,complib='zlib')	
        for di in sampler:
            df=pd.DataFrame(fp.tidy_trajectories(di))
            df['rep']=REPID*len(df.index)
            df.reset_index(['rep'],inplace=True,drop=True)
            output.append('trajectories',df)
            REPID+=1
        output.close()
    elif isinstance(sampler,fp.QtraitStatsSampler):
        ##Write in append mode
        output = pd.HDFStore(outputFilename,mode,complevel=6,complib='zlib')
        for s in sampler:
            df=pd.DataFrame(i)
            df['rep']=[REPID]*len(df.index)
            REPID+=1
        output.append('stats',pd.concat(data))
        output.close()
    else:
        raise RuntimeError("uh oh: sampler type not recognized for output.  We shouldn't have gotten this far!")
    return REPID
       
def make_parser():
    parser=argparse.ArgumentParser(description="Single region, simple demography, additive effects.")

    parser.add_argument('--seed','-S',type=int,default=0,metavar='SEED',help='Random number seed (default 0)')
    #parser.add_argument('--theta','-th',type=float,default=100.,metavar='THETA',help="4Nu/region (default 100.0)")
    parser.add_argument('--recrate','-r',type=float,default=0.5,help="Genetic distance of region.  Default = 0.5 = 50cM")
    parser.add_argument('--tsample','-t',type=int,default=0,metavar='TSAMPLE',help="How often to apply sampler. Default = 0 = never")
    parser.add_argument('--sigmu',type=float,default=0.25,metavar='SIGMU',help="Std. dev. of effect sizes. Default = 0.25")
    parser.add_argument('--H2',type=float,default=1.0,help="Desired heritability.  Default of 1.0. We use HoC approximation to determine std. dev. of environmental/random effects to give this heritability on average")
    parser.add_argument('--VS',type=float,default=1.0,help="V(S). Default to 1.0.")
    parser.add_argument('--mutrate','-m',type=float,default=-1.0,help="Mutation rate to selected variants. Required.",required=True)
    parser.add_argument("--outfile",'-O',type=str,default=None,help="Output file. Required",required=True)
    parser.add_argument('--trait',choices=TRAITS,default='additive',help="Genotype to phenotype model. Default is additive.")
    parser.add_argument('--sampler',choices=SAMPLERS,help="Temporal sampler type. Default is None")
    parser.add_argument('--dominance',type=float,default=1.0,help="Dominance. Default is 1.0 (co-dominance).")
    parser.add_argument('--ncores',type=int,default=64,help="No. simultaneous simulations to run. Default = 64")
    parser.add_argument('--nbatches',type=int,default=16,help="No. 'batches' of simulations to run.  Total number of replicates will be NCORES*NBATCHES")
    #Positional arguments
    parser.add_argument('popsize',type=int,default=10000)
    parser.add_argument('optimum',type=float,default=0.0)

    return parser

def main():
    parser=make_parser()
    parser.parse_args(sys.argv[1:])
    sys.exit(0) 
    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:e:H:S:O:N:t:s:F:r:n:d:",["theta=","rho=","trait=","sampler=","nsam=","cores=","batches="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    #set up default params
    N=1000   # pop size
    t=None     # 0.1N
    e = 0.25 # s.d. of effect sizes
    S = 1    # V(S)
    H = 1.0 # desired b-sense H^2
    m = None # Mutation rate (per gamete, per generation) to alleles affecting trait value
    r = 0.5 # rec. rate (per diploid, per gen)
    Opt = 0.0  # Value of optimum after 10N gens
    ofile=None
    seed = 0
    theta = 100 # mutational size of region where neutral mutations occur
    rho = 100   # recombination size of region where neutral mutations occur
    traitString="additive"
    samplerString=None
    ncores=64
    nbatches=16
    dominance=1.0
    for o,a in opts:
        if o == '-m':
            m = float(a)
        elif o == '-e':
            e = float(a)
        elif o == '-H':
            H = float(a)
        elif o == '-S':
            S = float(a)
        elif o == '-O':
            Opt = float(a)
        elif o == '-N':
            N=int(a)
        elif o == '-t':
            t = int(a)
        elif o == '-s':
            seed = int(a)
        elif o == '-F':
            ofile = a
        elif o == '-r':
            r = float(a)
        elif o == '--nsam':
            ssize = int(a)
        elif o == '--theta':
            theta=float(a)
        elif o == '--rho':
            rho=float(a)
        elif o == '--trait':
            traitString = a
            if a not in valid_trait_models:
                print "Error: invalid trait model"
                usage()
                sys.exit(0)
        elif o == '--sampler':
            samplerString = a
            if a not in valid_sampler_names:
                print "Error: invalid sampler name"
                usage()
                sys.exit(0)
        elif o == '--cores':
            ncores=int(a)
            if ncores < 1:
                print "--ncores must be > 0"
                usage()
                sys.exit(0)
        elif o == '--batches':
            nbatches=int(a)
            if nbatches < 1:
                print "--nbatches myst be > 0"
                usage()
                sys.exit(0)

    if samplerString is None:
        print "Error: sampler must be defined"
        usage()
        sys.exit(0)
        
    if t is None:
        print "Error: sampling interval must be defined"
        usage()
        sys.exit(0)
    if H is None:
        usage()
        sys.exit(2)
    if m is None:
        usage()
        sys.exit(2)
    if ofile is None:
        usage()
        sys.exit(2)

    ##Figure out sigma_E from params
    sigE = get_sigE_additive(m,S,H)

    trait = get_trait_model(traitString)

    nlist = np.array([N]*(10*N),dtype=np.uint32)
    rng=fp.GSLrng(seed)
    mu_neutral = 0.0
    nregions=[]
    sregions=[fp.GaussianS(0,1,1,e,dominance)]
    recregions=[fp.Region(0,1,1)]
    
    if samplerString == 'stats':
        output=pd.HDFStore(ofile,'w',complevel=6,complib='zlib')
        output.close()
    REPID=0
    for BATCH in range(nbatches):
        pops=fp.SpopVec(ncores,N)
        sampler=get_sampler_type(samplerString,len(pops),0.0)
        qt.evolve_regions_qtrait_sampler_fitness(rng,pops,sampler,trait,
                                                 nlist,
                                                 0.0,
                                                 m,
                                                 r,
                                                 nregions,
                                                 sregions,
                                                 recregions,
                                                 t,
                                                 sigE,
                                                 optimum=0.0,
                                                 VS=S)
        if samplerString == 'freq':
            dummy=write_output(sampler,ofile,REPID,BATCH,'w')
        elif samplerString == 'stats':
            dummy=write_output(sampler,ofile,REPID,BATCH,'a')
        sampler=get_sampler_type(samplerString,len(pops),Opt)
        qt.evolve_regions_qtrait_sampler_fitness(rng,pops,sampler,trait,
                                                 nlist,
                                                 0.0,
                                                 m,
                                                 r,
                                                 nregions,
                                                 sregions,
                                                 recregions,
                                                 t,
                                                 sigE,
                                                 optimum=Opt,
                                                 VS=S)
        if samplerString == 'freq':
            #Append this time!
            REPID=write_output(sampler,ofile,REPID,BATCH,'a')
        elif samplerString == 'stats':
            REPID=write_output(sampler,ofile,REPID,BATCH,'a')

if __name__ == "__main__":
    main()
