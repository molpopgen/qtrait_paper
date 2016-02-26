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
        opts, args = getopt.getopt(sys.argv[1:],"m:e:H:S:O:N:s:r:",["fixed=","lost="]))
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    #set up default params
    N=1000   # pop size
    e = 0.25 # s.d. of effect sizes
    S = 1    # V(S)
    H = None # desired b-sense H^2
    m = None # Mutation rate (per gamete, per generation) to alleles affecting trait value
    r = 0.5 # rec. rate (per diploid, per gen)
    Opt = 0.0  # Value of optimum after 10N gens
    fixationsFile=None
    lostFile=None
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
        elif o == '-s':
            seed = int(a)
        elif o == '-r':
            r = float(a)
        elif o == '--fixed':
            fixationsFile=a
        elif o == '--lost':
            lostFile=a

    if H is None:
        usage()
        sys.exit(2)
    if m is None:
        usage()
        sys.exit(2)
    if ofile is None:
        usage()
        sys.exit(2)
    if fixationsFile is None or lostFile is None:
        usage()
        sys.exit(2)

    rng = fp.GSLrng(seed)
    hdf_fixed = pd.HDFStore(fixationsFile,'w',complevel=6,complib='zlib')
    hdf_fixed.open()
    hdf_lost = pd.HDFStore(lostFile,'w',complevel=6,complib='zlib')
    hdf_list.open()
    
    sigE = get_sigE_additive(m,S,H)

    nregions = []
    recregions = [fp.Region(0,1,1)]
    sregions = [fp.GaussianS(0,1,1,e)]
    #population size over time -- constant & we re-use this over and over
    nlist = np.array([N]*(10*N),dtype=np.uint32)
    #16 batches of 64 runs = 1024 replicates  
    REPLICATE=0
    for i in range(16):
        #set up populations
        pop=fp.popvec(64,N)
        #Evolve to equilibrium
        traj1 = qt.evolve_qtrait_track(rng,
                                      pops,
                                      nlist[0:],
                                      0,
                                      m,
                                      r,
                                      nregions,sregions,recregions,
                                      sigmaE=sigE,
                                      track=1,
                                      VS=S) ##Do not track popstats during "burn-in"



        #evolve for another 3N generations at new optimum
        traj2=qt.evolve_qtrait_track(rng,pops,
                                     nlist[0:(3*N)],
                                     0,
                                     m,
                                     r,
                                     nregions,sregions,recregions,
                                     sigmaE=sigE,
                                     VS=S,optimum=Opt,track=1)
        AGES=[]
        FIXATIONS=[]
        for j in range(len(pop)):
            #Merge all trajectories for this replicate
            df = pd.concat([pd.DataFrame(traj1[j]),
                            pd.DataFrame(traj2[j]),
                            pd.DataFrame(traj3[j])])
            for name,group in df.groupby(['pos','esize']):
                if group.freq.max() < 1:  #mutation did not reach fixation...
                    if group.generation.max()-group.generation.min()>1: #... and it lived > 1 generation ...
                       AGES.append({'rep':REPLICATE,
                                    'esize':name[1],
                                    'origin':group.generation.min(),
                                    'final_g':group.generation.max(),
                                    'max_q':group.freq.max(),
                                    'last_q':group.freq.iloc[-1]})
                else: #mutation did reach fixation!
                    FIXATIONS.append({'rep':REPLICATE,
                                      'esize':name[1],
                                      'origin':group.generation.min(),
                                      'final_g':group.generation.max()})
            REPLICATE+=1

        hdf_fixed.append('fixations',pd.DataFrame(FIXATIONS))
        hdf_lost.append('allele_ages',pd.DataFrame(AGES))

    hdf_fixed.close()
    hdf_lost.close()

if __name__ == "__main__":
    main()
