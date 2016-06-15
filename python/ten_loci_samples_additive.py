#This script runs simulations and outputs various
#summaries of variability @ neutral markers

import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import libsequence.polytable as polyt
import KRT_qtrait_paper.summstatsParallel as sstats
import numpy as np
import pandas as pd
import getopt
import sys

def usage():
    print "something went wrong"

statNames=['tajd','thetaw','thetapi',"nd1",
           'hprime','nSL','iHS','H12','H1','H2H1']

def get_summstats_parallel(x):
    rv=[]
    for REP in x:
        for LOCUS in range(10):
            gen=[i[0] for i in REP]
            locus=[i[1][LOCUS]['genotypes'][0] for i in REP]
            v=sstats.simDataVec(locus)
            stats = sstats.getSummStatsParallel(v)
            DF=pd.DataFrame(stats)
            DF['gen']=gen
            rv.append(DF)
    return pd.concat(rv)
#    for GEN in range(len(x[0])):
#        for LOCUS in range(10):
#        for GEN in range
#            locus1=[i[0][1][LOCUS]['genotypes'][0] for i in x[GEN]]
#            v=sstats.simDataVec(locus1)
#            stats = sstats.getSummStatsParallel(v)
#            print stats
# def get_summstats(x):
#     gen=x[0] #This is the generation
#     rv=[]
#     LOCUS=0
#     for i in x[1]:
#         sd=polyt.simData(i['genotypes'][0])
#         ps=sstats.polySIM(sd)
#         nSL=sstats.std_nSLiHS(sd)
#         g=sstats.garudStats(sd)
#         values=[ps.tajimasd(),
#                 ps.thetaw(),
#                 ps.thetapi(),
#                 ps.numexternalmutations(),
#                 ps.hprime(),nSL[0],nSL[1],
#                 g['H12'],g['H1'],g['H2H1']]
#         rv.append(pd.DataFrame({'gen':[gen]*len(values),'locus':[LOCUS]*len(values),'stats':statNames,'values':values}))
#         LOCUS += 1
#     return pd.concat(rv)

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
    nlist=np.array([N]*(100),dtype=np.uint32)
    for BATCH in range(16): #16*64=1024
        x = fp.popvec_mloc(NREPS,N,NLOCI)

        samples = qtm.evolve_qtraits_mloc_sample(rnge,rngs,x,nlist,
                                                 [mu_n_region]*NLOCI,
                                                 [mu_del_ttl/float(NLOCI)]*NLOCI,
                                                 [e]*NLOCI,
                                                 [little_r_per_locus]*NLOCI,
                                                 [0.5]*(NLOCI-1),#loci unlinked
                                                 sample=t,nsam=nsam,VS=S)

        RTMP=REP
        print "done"
        print get_summstats_parallel(samples)
        #for si in samples:
        #    ti=pd.concat([get_summstats(i) for i in si if i[0] != 10*N])
        #    out.append('summstats',ti)
        #    RTMP+=1

        samples = qtm.evolve_qtraits_mloc_sample(rnge,rngs,x,nlist,
                                                 [mu_n_region]*NLOCI,
                                                 [mu_del_ttl/float(NLOCI)]*NLOCI,
                                                 [e]*NLOCI,
                                                 [little_r_per_locus]*NLOCI,
                                                 [0.5]*(NLOCI-1),#loci unlinked
                                                 sample=t,nsam=nsam,VS=S,optimum=Opt)
        
        #for si in samples:
        #    ti=pd.concat([get_summstats(i) for i in si])
        #    out.append('summstats',ti)
        #    REP+=1

    out.close()

if __name__ == "__main__":
    main()
