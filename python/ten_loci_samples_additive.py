#This script runs simulations and outputs various
#summaries of variability @ neutral markers

#Note: this script depends on KT_qtrait_paper.summstatsParallel,
#which is a Cython extension module extending pylibseq.

import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import libsequence.polytable as polyt
import libsequence.extensions as ext
import KRT_qtrait_paper.summstatsParallel as sstats
import numpy as np
import pandas as pd
import getopt
import sys

def usage():
    print "something went wrong"

statNames=['tajd','thetaw','thetapi',"nd1",
           'hprime','nSL','iHS','H12','H1','H2H1']

def get_summstats_parallel(x,REPID,out):
    rv=[]
    for REP in x:
        for LOCUS in range(10):
            gen=[i[0] for i in REP]
            locus=[i[1][LOCUS]['genotypes'][0] for i in REP]
            v=ext.SimDataVec(locus)
            stats = sstats.getSummStatsParallel(v)
            DF=[pd.DataFrame(i.items(),columns=['stat','value']) for i in stats]
            for i in range(len(DF)):
                DF[i]['gen']=[gen[i]]*len(DF[i].index)
                DF[i]['locus']=[LOCUS]*len(DF[i].index)
                DF[i]['replicate']=[REPID]*len(DF[i].index)
            rv.append(pd.concat(DF))
        REPID+=1
    out.append('summstats',pd.concat(rv))

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:e:H:S:O:N:t:s:r:F:n:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    N=1000   # pop size
    e = 0.25 # s.d. of effect sizes
    S = 1    # V(S)
    H = None # desired b-sense H^2
    mu_del_ttl = None # Mutation rate (per gamete, per generation) to alleles affecting trait value
    r = 0.5 # rec. rate b/w loci (per diploid, per gen)
    Opt = 0.0  # Value of optimum after 10N gens
    ofile=None
    seed = 0
    t=None
    nsam=100
    for o,a in opts:
        if o == '-m':
            mu_del_ttl = float(a)
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
        elif o == '-s':
            seed = int(a)
        elif o == '-r':
            r = float(a)
        elif o == '-F':
            ofile=a
        elif o == '-t':
            t=int(a)
        elif o == '-n':
            nsam=int(a)

    if H is None:
        usage()
        sys.exit(2)
    if mu_del_ttl is None:
        usage()
        sys.exit(2)
    if ofile is None:
        usage()
        sys.exit(2)
    if t is None:
        t=0.1*float(N)
    #Constants:
    NLOCI=10
    NREPS=64
    ##The next 2 are the 'sizes' of each locus in scaled params.  Roughly 100kb in 'humans' in terms of pi.
    theta=100.0
    rho=100.0

    #Can start working now:
    REP=0
    out=pd.HDFStore(ofile,"w",complevel=6,complib='zlib')

    little_r_per_locus = rho/(4.0*float(N))
    mu_n_region=theta/(4.0*float(N))
    
    rnge=fp.GSLrng(seed)
    rngs=fp.GSLrng(seed)
    nlist=np.array([N]*(10*N),dtype=np.uint32)
    fitness = qtm.MlocusAdditiveTrait()
    sregions=[fp.GaussianS(0,1,1,e,1.0)]*NLOCI
    for BATCH in range(16): #16*64=1024
        x = fp.MlocusPopVec(NREPS,N,NLOCI)
        sampler = fp.PopSampler(len(x),nsam,rngs)
        qtm.evolve_qtraits_mloc_sample_fitness(rnge,x,sampler,fitness,nlist,
                                               [mu_n_region]*NLOCI,
                                               [mu_del_ttl/float(NLOCI)]*NLOCI,
                                               sregions,
                                               [little_r_per_locus]*NLOCI,
                                               [0.5]*(NLOCI-1),#loci unlinked
                                               sample=t,VS=S)
        samples=sampler.get()
        get_summstats_parallel(samples,REP,out)
        sampler = fp.PopSampler(nsam,rngs)
        qtm.evolve_qtraits_mloc_sample_fitnes(rnge,x,sampler,fitness,nlist,
                                              [mu_n_region]*NLOCI,
                                              [mu_del_ttl/float(NLOCI)]*NLOCI,
                                              sregions,
                                              [little_r_per_locus]*NLOCI,
                                              [0.5]*(NLOCI-1),#loci unlinked
                                              sample=t,VS=S,optimum=Opt)
        samples=sampler.get()
        get_summstats_parallel(samples,REP,out)
        REP += NREPS
    out.close()

if __name__ == "__main__":
    main()
