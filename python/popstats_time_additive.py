import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import pandas as pd,getopt,sys,warnings,math

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

    ##Figure out sigma_E from params
    EVG = 4.0*m*S*S
    EVE = EVG*(1.0-H)/H
    sigE = math.sqrt(EVE)

    rng = fp.GSLrng(seed)
    hdf = pd.HDFStore(ofile,'w',complevel=6,complib='zlib')
    hdf.open()
    #16 batches of 64 runs = 1024 replicates
    for i in range(16):
        #Simulate 10N-1 generations to equilbrium with optimum of 0
        nlist = np.array([N]*(10*N - 1),dtype=np.uint32)
        #Evolve to equilibrium
        pops = qt.evolve_qtrait(rng,
                                64,
                                N,
                                nlist[0:],
                                0,
                                m,
                                r,
                                [],
                                [fp.GaussianS(0,1,1,e)],
                                [fp.Region(0,1,1)],
                                sigmaE=sigE,
                                VS=S) ##Do not track popstats during "burn-in"

        #simulate another 2*N generations, sampling stats every 't' generations
        nlist = np.array([N]*(2*N),dtype=np.uint32)

        stats = qt.evolve_qtrait_popstats(rng,pops,
                                          nlist[0:],
                                          0,
                                          m,
                                          r,
                                          [],
                                          [fp.GaussianS(0,1,1,e)],
                                          [fp.Region(0,1,1)],
                                          sigmaE=sigE,
                                          trackStats=t,
                                          VS=S)
        for j in range(len(pops)):
            hdf.append('popstats',pd.DataFrame(stats[j]))

        #Now, shift the optimum, and evolve for another 3N generations,
        #sampling every t generations

        #We shift the optimim and sample immediately, so that we get
        #any "spikes" in VG, etc.
        nlist = np.array([N]*1,dtype=np.uint32)
        stats=qt.evolve_qtrait_popstats(rng,pops,
                              nlist[0:],
                              0,
                              m,
                              r,
                              [],
                              [fp.GaussianS(0,1,e,1)],
                              [fp.Region(0,1,1)],
                              sigmaE=sigE,
                              VS=S,optimum=Opt,trackStats=t)

        for j in range(len(pops)):
            hdf.append('popstats',pd.DataFrame(stats[j]))

        #evolve for another 3N generations post optimum-shift
        nlist = np.array([N]*(3*N),dtype=np.uint32)

        stats=qt.evolve_qtrait_popstats(rng,pops,
                                  nlist[0:],
                                  0,
                                  m,
                                  r,
                                  [],
                                  [fp.GaussianS(0,1,1,e)],
                                  [fp.Region(0,1,1)],
                                  sigmaE=sigE,
                                  VS=S,optimum=Opt,trackStats=t)

        for j in range(len(pops)):
            hdf.append('popstats',pd.DataFrame(stats[j]))

    hdf.close()

if __name__ == "__main__":
    main()
