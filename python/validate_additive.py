import fwdpy as fp
import fwdpy.qtrait as qt
import pandas as pd
import numpy as np

hdf = pd.HDFStore("popstats.h5",'w')
hdf.open()
k=0

N=1000
nlist = np.array([N]*10*N,dtype=np.uint32)

nregions=[]
##Gassian with sd=0.1
recregions=[fp.Region(0,1,1)]

muvals=[1e-4,1e-3,2.5e-3,5e-3,7.5e-3,1e-2]
sigmuvals=[0.125,0.25,0.5]
sigevals=[0.1,0.25]
recvals=[0.,1./(4.*float(N)),10./(4.*float(N)),100./(4.*float(N))]

rng = fp.GSLrng(27288)

BATCHES=8
for mu in muvals:
    for sigmu in sigmuvals:
        sregions = [fp.GaussianS(0,1,1,sigmu)]
        for sige in sigevals:
            eh2 = 4.*mu/(4*mu + sige**sige)
            for rec in recvals:
                for i in range(BATCHES):
                    pops = qt.evolve_qtrait(rng,64,
                                            N,
                                            nlist[0:],
                                            0,
                                            mu,
                                            rec,
                                            nregions,
                                            sregions,
                                            recregions,
                                            sige)
                    stats = fp.diploid_traits(pops)
                    vals=[]
                    for s in stats:
                        G = np.array([si['g'] for si in s],dtype=np.float64)
                        val = {'rep':k,
                               'mu':mu,
                               'sige':sige,
                               'r':rec,
                               'sigmu':sigmu,
                               'VG':np.var(G)}
                        vals.append(val)
                        k=k+1

                    hdf.append('popstats',pd.DataFrame(vals))
hdf.close()
                        
