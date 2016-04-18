import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import pandas as pd,getopt,sys,warnings,math

##import local functions
from single_region_common import *

#corresponds to making Figure 4 of the grant, which is VG, etc. over time

def usage():
    print "something went wrong"

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:e:H:S:O:N:t:s:F:r:")
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
    REPLICATE=0 
    for i in range(16):
        pops=fp.popvec(64,N)
        #Evolve to equilibrium, tracking along the way
        stats = qt.evolve_qtrait_popstats(rng,
                                         pops,
                                         nlist[0:],
                                         0,
                                         m,
                                         r,
                                         nregions,sregions,recregions,
                                         sigmaE=sigE,
                                         trackStats=t,
                                         VS=S) 
        RTEMP=REPLICATE
        for si in stats:
            t=pd.DataFrame(si)
            t['rep']=[RTEMP]*len(t.index)
            RTEMP+=1
            hdf.append('popstats',t)
        #simulate another 10*N generations, sampling stats every 't' generations
        stats = qt.evolve_qtrait_popstats(rng,pops,
                                          nlist[0:],
                                          0,
                                          m,
                                          r,
                                          nregions,sregions,recregions,
                                          sigmaE=sigE,
                                          trackStats=t,
                                          VS=S,optimum=Opt)
        for si in stats:
            t=pd.DataFrame(si)
            t['rep']=[REPLICATE]*len(t.index)
            hdf.append('popstats',t)
            REPLICATE+=1
        #for j in range(len(pops)):
        #    hdf.append('popstats',pd.DataFrame(stats[j]))

    hdf.close()

if __name__ == "__main__":
    main()
