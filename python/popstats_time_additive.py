import fwdpy as fp
import fwdpy.qtrait as qt
#import fwdpy.fitness as fpw
import numpy as np
import pandas as pd,getopt,sys,warnings,math

##import local functions
from single_region_common import *

#corresponds to making Figure 4 of the grant, which is VG, etc. over time

def usage():
    print "something went wrong"

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:e:H:S:O:N:t:s:F:r:",["cores=","batches="])
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
    ncores=64
    nbatches=16
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
        elif o == '--cores':
            ncores=int(a)
        elif o == '--batches':
            nbatches=int(a)
            
    if t is None:
        t = int(0.1*float(N))
    if H is None:
        usage()
        sys.exit(2)
    if m is None:
        usage()
        sys.exit(2)
    if ofile is None:
        usage()
        sys.exit(2)



    rng = fp.GSLrng(seed)
    hdf = pd.HDFStore(ofile,'w',complevel=6,complib='zlib')
    hdf.open()
    
    sigE = get_sigE_additive(m,S,H)

    nregions = []
    recregions = [fp.Region(0,1,1)]
    sregions = [fp.GaussianS(0,1,1,e)]
    #population size over time -- constant & we re-use this over and over
    nlist = np.array([N]*(10*N),dtype=np.uint32)
    #16 batches of 64 runs = 1024 replicates
    fitness=qt.SpopAdditiveTrait()
    REPLICATE=0 
    for i in range(nbatches):
        pops=fp.SpopVec(ncores,N)
        sampler = fp.QtraitStatsSampler(len(pops),0.0)
        #Evolve to equilibrium, tracking along the way
        qt.evolve_regions_qtrait_sampler_fitness(rng,
                                                 pops,
                                                 sampler,fitness,
                                                 nlist[0:],
                                                 0,
                                                 m,
                                                 r,
                                                 nregions,sregions,recregions,
                                                 sigmaE=sigE,
                                                 sample=t,
                                                 VS=S,optimum=0)
        stats=sampler.get()
        RTEMP=REPLICATE
        for si in stats:
            ti=pd.DataFrame(si)
            ti['rep']=[RTEMP]*len(ti.index)
            RTEMP+=1
            hdf.append('popstats',ti)
        #simulate another 10*N generations, sampling stats every 't' generations
        sampler = fp.QtraitStatsSampler(len(pops),Opt)
        qt.evolve_regions_qtrait_sampler_fitness(rng,
                                                 pops,
                                                 sampler,fitness,
                                                 nlist[0:],
                                                 0,
                                                 m,
                                                 r,
                                                 nregions,sregions,recregions,
                                                 sigmaE=sigE,
                                                 sample=t,
                                                 VS=S,optimum=Opt)
        stats=sampler.get()
        for si in stats:
            ti=pd.DataFrame(si)
            ti['rep']=[REPLICATE]*len(ti.index)
            hdf.append('popstats',ti)
            REPLICATE+=1

    hdf.close()

if __name__ == "__main__":
    main()
