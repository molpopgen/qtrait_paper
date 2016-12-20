import fwdpy as fp
import fwdpy.qtrait as qt
import pandas as pd
import numpy as np


hdf = pd.HDFStore("validation_popstats.h5",'w',complevel=6,complib='zlib')
hdf2 = pd.HDFStore("validation_popstats.h5",'w',complevel=6,complib='zlib')

N=5000
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
            eh2 = 4.*mu/(4.*mu + sige**2.)
            for rec in recvals:
		k=0
                for i in range(BATCHES):
                    pops = qt.evolve_reqions_qtrait(rng,64,
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
                    views = fp.view_mutations(pops)
                    hdf.append('popstats',pd.DataFrame(vals))
hdf.close()
                        
