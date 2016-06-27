import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import gzip,getopt,sys,warnings,math
import pandas as pd
##import local functions
from single_region_common import *

valid_trait_models=['additive','mult']
trait_models = {'additive':qt.SpopAdditiveTrait(),
                'mult':qt.SpopMultTrait()}

valid_sampler_names=['stats','freq']

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

def write_output(sampler,output,REPID):
    if isinstance(sampler,fp.FreqSampler):
        data=[pd.DataFrame(i) for i in fp.tidy_trajectories(sampler.get())]
        for df in data:
            df['rep']=[REPID]*len(df.index)
            REPID+=1
        output.append('trajectories',pd.concat(data))
    elif isinstance(sampler,fp.QtraitStatsSampler):
        data=[pd.DataFrame(i) for i in sampler.get()]
        for df in data:
            df['rep']=[REPID]*len(df.index)
            REPID+=1
        output.append('stats',pd.concat(data))
    else:
        raise RuntimeError("uh oh: sampler type not recognized for output.  We shouldn't have gotten this far!")
    return REPID
        
def main():
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
    H = None # desired b-sense H^2
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
    output=pd.HDFStore(ofile,'w',complevel=6,complib='zlib')
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
        ff=pd.HDFStore("first.h5",'w')
        X=fp.tidy_trajectories(sampler.get())
        ff.append('traj',pd.DataFrame(X[0]))
        ff.close()
        dummy=write_output(sampler,output,REPID)
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
        ff=pd.HDFStore("last.h5",'w')
        X=fp.tidy_trajectories(sampler.get())
        ff.append('traj',pd.DataFrame(X[0]))
        ff.close()
        REPID=write_output(sampler,output,REPID)

    output.close()
if __name__ == "__main__":
    main()
