import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import numpy as np
import pandas as pd
import getopt
import sys

def usage():
    print "something went wrong"


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:e:H:S:O:N:t:s:r:F:")
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
        t=int(0.1*float(N))
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
    nlist=np.array([N]*(10*N),dtype=np.uint32)
    fitness = qtm.MlocusAdditiveTrait()
    sregions=[fp.GaussianS(0,1,1,e,1.0)]*NLOCI
    for BATCH in range(16): #16*64=1024
        x = fp.MlocusPopVec(NREPS,N,NLOCI)
        sampler = fp.QtraitStatsSampler(len(x),0.0)
        qtm.evolve_qtraits_mloc_sample_fitness(rnge,x,sampler,fitness,nlist,
                                               [mu_n_region]*NLOCI,
                                               [mu_del_ttl/float(NLOCI)]*NLOCI,
                                               sregions,
                                               [little_r_per_locus]*NLOCI,
                                               [0.5]*(NLOCI-1),#loci unlinked
                                               sample=t,VS=S)

        RTMP=REP
        stats=sampler.get()
        for si in stats:
            tdf=pd.DataFrame([i for i in si if i['generation'] != 10*N])
            tdf['rep']=[RTMP]*len(tdf.index)
            out.append('popstats',tdf)
            RTMP+=1

        sampler = fp.QtraitStatsSampler(len(x),Opt)
        qtm.evolve_qtraits_mloc_sample_fitness(rnge,x,sampler,fitness,nlist,
                                               [mu_n_region]*NLOCI,
                                               [mu_del_ttl/float(NLOCI)]*NLOCI,
                                               sregions,
                                               [little_r_per_locus]*NLOCI,
                                               [0.5]*(NLOCI-1),#loci unlinked
                                               sample=t,VS=S,optimum=Opt)
        stats=sampler.get()
        for si in stats:
            tdf=pd.DataFrame(si)
            tdf['rep']=[REP]*len(tdf.index)
            out.append('popstats',tdf)
            REP+=1

    out.close()

if __name__ == "__main__":
    main()
