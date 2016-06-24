import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import numpy as np
import pandas as pd
#import matplotlib
#import matplotlib.pyplot as plt
import copy,sys

N=1000
NLOCI=10
NREPS=40

rnge=fp.GSLrng(100)
nlist=np.array([N]*(15*N),dtype=np.uint32)
theta_neutral_per_locus=0.0
rho_per_locus=100.0
little_r_per_locus=rho_per_locus/(4.0*float(N))
mu_n_region=theta_neutral_per_locus/(4.0*float(N))
mu_del_ttl=1e-3

sigmas=[0.1]*NLOCI

hdf=pd.HDFStore("cumVA.h5",'w')

REPLICATE=0
for batch in range(125): #4,000 reps...
    x = qtm.evolve_qtraits_mloc(rnge,NREPS,N,NLOCI,
                                nlist[:(10*N)],
                                [mu_n_region]*NLOCI,
                                [mu_del_ttl/float(NLOCI)]*NLOCI,
                                sigmas,
                                [little_r_per_locus]*NLOCI,
                                [0.5]*(NLOCI-1),#loci unlinked
                                )
    traj2=qtm.evolve_qtraits_mloc_VA(rnge,x,nlist,  #Only evolve N further generations...
                                     [mu_n_region]*NLOCI,
                                     [mu_del_ttl/float(NLOCI)]*NLOCI,
                                     sigmas,
                                     [little_r_per_locus]*NLOCI,
                                     [0.5]*(NLOCI-1),#loci unlinked
                                     sample=10,optimum=0.5)
    
    # df=[pd.DataFrame(i) for i in traj]
    # TEMP=REPLICATE
    # for i in df:
    #     i['rep']=[TEMP]*len(i.index)
    #     TEMP+=1
    # df=pd.concat(df)
    # df=df[df.generation<10000] ##Remove time of optimum shift, so that it only happens once in the output
    # hdf.append('cumVA',df)
    df=[pd.DataFrame(i) for i in traj2]
    TEMP=REPLICATE
    for i in df:
        i['rep']=[TEMP]*len(i.index)
        TEMP+=1
    hdf.append('cumVA',pd.concat(df))
    REPLICATE+=NREPS

