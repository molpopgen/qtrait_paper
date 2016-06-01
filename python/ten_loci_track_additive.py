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
        opts, args = getopt.getopt(sys.argv[1:],"m:e:H:S:O:N:s:r:",["fixed=","ages=","traj="])
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
    fixationsFile=None
    lostFile=None
    trajFile=None
    seed = 0
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
        elif o == '--fixed':
            fixationsFile=a
        elif o == '--ages':
            lostFile=a
        elif o == '--traj':
            trajFile=a

    if H is None:
        usage()
        sys.exit(2)
    if mu_del_ttl is None:
        usage()
        sys.exit(2)
    if trajFile is None:
        usage()
        sys.exit(2)
        
    #Constants:
    NLOCI=10
    NREPS=64
    theta=100.0
    rho=100.0

    #Can start working now:
    REP=0
    out=pd.HDFStore(trajFile,"w")

    little_r_per_locus = rho/(4.0*float(N))
    mu_n_region=theta/(4.0*float(N))
    
    rnge=fp.GSLrng(seed)
    nlist=np.array([N]*(10*N),dtype=np.uint32)
    for BATCH in range(16): #16*64=1024
        x = fp.popvec_mloc(NREPS,N,NLOCI)

        traj1 = qtm.evolve_qtraits_mloc_track(rnge,x,nlist,
                                              [mu_n_region]*NLOCI,
                                              [mu_del_ttl/float(NLOCI)]*NLOCI,
                                              [e]*NLOCI,
                                              [little_r_per_locus]*NLOCI,
                                              [0.5]*(NLOCI-1),#loci unlinked
                                              1,VS=S)

    
        traj2 = qtm.evolve_qtraits_mloc_track(rnge,x,nlist,
                                              [mu_n_region]*NLOCI,
                                              [mu_del_ttl/float(NLOCI)]*NLOCI,
                                              [e]*NLOCI,
                                              [little_r_per_locus]*NLOCI,
                                              [0.5]*(NLOCI-1),#loci unlinked
                                              1,VS=S,optimum=Opt)
        
        m=fp.tidy_trajectories(fp.merge_trajectories(traj1,traj2))
        for i in m:
            temp=pd.DataFrame(i)
            temp['rep']=[REP]*len(temp.index)
            out.append('traj',temp)
            REP+=1

    out.close()

if __name__ == "__main__":
    main()
