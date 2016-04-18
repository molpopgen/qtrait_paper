import fwdpy as fp
import fwdpy.qtrait as qt
import numpy as np
import pandas as pd,getopt,sys,warnings,math
import copy
import feather
##import local functions
from single_region_common import *

#corresponds to making Figure 4 of the grant, which is VG, etc. over time

def usage():
    print "something went wrong"

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:e:H:S:O:N:s:r:",["fixed=","ages=","traj="])
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
    trajFile=None
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
        elif o == '--ages':
            lostFile=a
        elif o == '--traj':
            trajFile=a

    if H is None:
        usage()
        sys.exit(2)
    if m is None:
        usage()
        sys.exit(2)
    if fixationsFile is None or lostFile is None or trajFile is None:
        usage()
        sys.exit(2)

    rng = fp.GSLrng(seed)
    hdf_fixed = pd.HDFStore(fixationsFile,'w',complevel=6,complib='zlib')
    hdf_fixed.open()
    hdf_lost = pd.HDFStore(lostFile,'w',complevel=6,complib='zlib')
    hdf_lost.open()
    hdf_traj = pd.HDFStore(trajFile,'w',complevel=6,complib='zlib')
    hdf_traj.open()
    sigE = get_sigE_additive(m,S,H)

    nregions = []
    recregions = [fp.Region(0,1,1)]
    sregions = [fp.GaussianS(0,1,1,e)]
    #population size over time -- constant & we re-use this over and over
    nlist = np.array([N]*(10*N),dtype=np.uint32)

    REPLICATE=0
    #16 batches of 64 runs = 1024 replicates
    for i in range(16):
        #set up populations
        pops=fp.popvec(64,N)
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

        #evolve for another 10N generations at new optimum
        traj2=qt.evolve_qtrait_track(rng,pops,
                                     nlist[0:],
                                     0,
                                     m,
                                     r,
                                     nregions,sregions,recregions,
                                     sigmaE=sigE,
                                     VS=S,optimum=Opt,track=1)
        AGES=[]
        FIXATIONS=[]
        #merge trajectories and get allele ages (parallelized via open MP)
        traj1=fp.merge_trajectories(traj1,traj2)

        ages = fp.allele_ages(traj1)
        REPTEMP=REPLICATE
        for ai in range(len(ages)):
            dfi=pd.DataFrame(ages[ai])
            dfi['rep']=[REPTEMP]*len(dfi.index)
            FIXATIONS.append(dfi[dfi['max_freq']==1.0])
            AGES.append(dfi[dfi['max_freq']<1.0])
            REPTEMP+=1
        # for j in range(len(pops)):
        #     #Merge all trajectories for this replicate
        #     df = pd.concat([pd.DataFrame(traj1[j]),
        #                     pd.DataFrame(traj2[j])])
        #     for name,group in df.groupby(['pos','esize']):
        #         if group.freq.max() < 1:  #mutation did not reach fixation...
        #             if group.generation.max()-group.generation.min()>1: #... and it lived > 1 generation ...
        #                AGES.append({'rep':REPLICATE,
        #                             'esize':name[1],
        #                             'origin':group.generation.min(),
        #                             'final_g':group.generation.max(),
        #                             'max_q':group.freq.max(),
        #                             'last_q':group.freq.iloc[-1]})
        #         else: #mutation did reach fixation!
        #             FIXATIONS.append({'rep':REPLICATE,
        #                               'esize':name[1],
        #                               'origin':group.generation.min(),
        #                               'final_g':group.generation.max()})
        #     REPLICATE+=1

        ##Add more info into the trajectories
        for t in traj1:
            LD=[]
            for i in t:
                I=int(0)
                for j in i[1]:
                    x=copy.deepcopy(i[0])
                    x['freq']=j
                    x['generation']=i[0]['origin']+I
                    I+=1
                    LD.append(x)
            d=pd.DataFrame(LD)
            d['rep']=[REPLICATE]*len(d.index)
            REPLICATE+=1
            hdf_traj.append('trajectories',d)

        hdf_fixed.append('fixations',pd.concat(FIXATIONS))
        hdf_lost.append('allele_ages',pd.concat(AGES))

    hdf_fixed.close()
    hdf_lost.close()
    hdf_traj.close()

if __name__ == "__main__":
    main()
