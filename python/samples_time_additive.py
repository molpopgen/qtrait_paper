#Track samples of size n at regular intervals,
#and write to a "pickled" archive

import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
try:
    import cPickle as pickle
except ImportError:
    import pickle

import gzip,getopt,sys,warnings,math

##import local functions
from single_region_common import *


def usage():
    print "something went wrong"

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:e:H:S:O:N:t:s:F:r:n:",["theta=","rho="])
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
    ssize = 50 #sample size
    theta = 100 # mutational size of region where neutral mutations occur
    rho = 100   # recombination size of region where neutral mutations occur
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
        elif o == '-n':
            ssize = int(a)
        elif o == '--theta':
            theta=float(a)
        elif o == '--rho':
            rho=float(a)

    if t is None:
        t = 0
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

    rng = fp.GSLrng(seed)

    PIK=gzip.open(ofile,"wb")
    #16 batches of 64 runs = 1024 replicates
    
    mun = theta/float(4*N)
    littler = rho/float(4*N)

    nreg=make_neutral_region()
    recstuff=make_buried_rec_region(littler)
    sreg=make_Gaussian_sregion(recstuff['beg'],recstuff['end'],sigma=e)
    rrg=recstuff['region']

    nlist = np.array([N]*(10*N),dtype=np.uint32)
    for i in range(16):
        #Evolve to equilibrium
        pops = qt.evolve_qtrait(rng,
                                64,
                                N,
                                nlist[0:((10*N)-1)],
                                mun,
                                m,
                                r,
                                nreg,sreg,rrg,
                                sigmaE=sigE,
                                VS=S) ##Do not track popstats during "burn-in"

        #simulate another 2*N generations, sampling stats every 't' generations
        samples = qt.evolve_qtrait_sample(rng,pops,
                                          nlist[0:(2*N)],
                                          mun,
                                          m,
                                          r,
                                          nreg,sreg,rrg,
                                          sigmaE=sigE,
                                          trackSamples=t,nsam=ssize,
                                          VS=S)
        for j in samples:
            pickle.dump(j,PIK)

        #We shift the optimum and sample immediately, so that we get
        #any "spikes" in VG, etc.
        samples =qt.evolve_qtrait_sample(rng,pops,
                                         nlist[0:1],
                                         0,
                                         m,
                                         r,
                                         [],
                                         [fp.GaussianS(0,1,e,1)],
                                         [fp.Region(0,1,1)],
                                         sigmaE=sigE,
                                         VS=S,optimum=Opt,trackSamples=t,nsam=ssize)
        for j in samples:
            pickle.dump(j,PIK)
 
        #evolve for another 3N generations post optimum-shift
        samples=qt.evolve_qtrait_sample(rng,pops,
                                        nlist[0:(3*N)],
                                        mun,
                                        m,
                                        r,
                                        nreg,sreg,rrg,
                                        sigmaE=sigE,
                                        VS=S,optimum=Opt,trackSamples=t,nsam=ssize)
        for j in samples:
            pickle.dump(j,PIK)

    PIK.close()

if __name__ == "__main__":
    main()
