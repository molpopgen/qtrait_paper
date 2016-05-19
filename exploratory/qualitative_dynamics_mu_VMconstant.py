#Generate a series of plots of VG, etc.,
#over time, to show relation of simulation
#params to HoC approximation

#This script varies mu with VS=1 and varies sigma_mu
#such that mu*(sigma_mu^2) is constant

import fwdpy as fp
import fwdpy.qtrait as qt
import pandas as pd
import numpy as np
import math,sys
import matplotlib
import matplotlib.pyplot as plt

nregions=[]
rregions=[fp.Region(0,1,1)]

reference_mu = 1e-3
recrate = 0.5
N=1000
simlen=10*N
Nlist = np.array([N]*(simlen),dtype=np.uint32)
#grid of VS values.  1 is our reference value
relative_mu=[0.25,0.5,1.0,5.0,10.0]

reference_sigma=math.sqrt(0.05)
reference_vm = 2.0*reference_mu*(math.pow(reference_sigma,2.0))
reference_vg = 4*reference_mu*1.0
rng=fp.GSLrng(1525152)

#plot VG, ebar, tbar,max_expl over time
fig=plt.figure(figsize=(10,10))
tbar=plt.subplot2grid((3,2),(1,0),rowspan=1,colspan=2)
VG=plt.subplot2grid((3,2),(0,0),rowspan=1,colspan=2)
sfs=plt.subplot2grid((3,2),(2,0),rowspan=1,colspan=2)

#hdf5 files to save the output for further plotting...
popstats_output = pd.HDFStore('popstats_vary_mu_VMconstant.h5','w',complevel=6,complib='zlib')
sfs_output = pd.HDFStore('sfs_vary_mu_VMconstant.h5','w',complevel=6,complib='zlib')

for rmu in relative_mu:
    mu = rmu*reference_mu
    sigma = math.sqrt(reference_vm/(2.0*mu))
    sregions=[fp.GaussianS(0,1,1,sigma)]
    vs=1.0
    #Generate 1e3 replicates with these params
    replabel=0
    DFS=[]
    sfs_at_opt=pd.Series([0.0]*(2*N))
    for batch in range(25):
        pops=fp.popvec(40,N)
        stats=qt.evolve_qtrait_popstats(rng,pops,
                                        Nlist[0:],
                                        0,
                                        mu,
                                        recrate,
                                        nregions,
                                        sregions,
                                        rregions,
                                        0.0, #No environmental noise
                                        trackStats=1,
                                        optimum=0.,
                                        VS=vs)
        muts=[fp.view_mutations(i) for i in pops]
        for mi in muts:
            for i in mi:
                sfs_at_opt[i['n']-1]+=1.0

        df=[pd.DataFrame([i for i in j if i['generation']!=simlen]) for j in stats]
        REP=replabel
        for i in df:
            i['rep']=[REP]*len(i.index)
            REP+=1
        DFS.extend(df)

        stats=qt.evolve_qtrait_popstats(rng,pops,
                                        Nlist[0:],
                                        0,
                                        mu,
                                        recrate,
                                        nregions,
                                        sregions,
                                        rregions,
                                        0.0, #No environmental noise
                                        trackStats=1,
                                        optimum=0.5*math.sqrt(vs),
                                        VS=vs)
        df=[pd.DataFrame(j) for j in stats]
        REP=replabel
        for i in df:
            i['rep']=[REP]*len(i.index)
            REP+=1
        DFS.extend(df)

        replabel+=40

    sfs_at_opt=sfs_at_opt.divide(float(replabel))
    df=pd.concat(DFS)

    dfg=df.groupby(['generation','stat']).mean()
    dfg.reset_index(inplace=True)
    dfg['mu']=[mu]*len(dfg.index)
    dfg['sigmu']=[sigma]*len(dfg.index)
    dfg['VS']=[vs]*len(dfg.index)

    #Write the means out to file
    popstats_output.append('popstats',dfg)
    sfsDF=pd.DataFrame()
    sfsDF['mu']=[mu]*len(sfs_at_opt)
    sfsDF['sigmu']=[sigma]*len(sfs_at_opt)
    sfsDF['VS']=[vs]*len(sfs_at_opt)
    sfsDF['msfs']=sfs_at_opt
    sfs_output.append('sfs',sfsDF)
    VG.plot(dfg[dfg.stat=='VG'].generation,dfg[dfg.stat=='VG'].value,
            label=r'$VS = $'+'{0:0.3}'.format(vs)+r', $\mu = $'+'{0:0.2}'.format(mu)+r', $\sigma_\mu = $'+'{0:0.2}'.format(sigma))
    tbar.plot(dfg[dfg.stat=='tbar'].generation,dfg[dfg.stat=='tbar'].value.divide(0.5*math.sqrt(vs)),
              label=r'$VS = $'+'{0:0.3}'.format(vs)+r', $\mu = $'+'{0:0.2}'.format(mu)+r', $\sigma_\mu = $'+'{0:0.2}'.format(sigma))
    sfs_at_opt=sfs_at_opt.divide(sfs_at_opt.sum())
    sfs.plot(np.arange(1.,11.,1.),sfs_at_opt[:10],
             label=r'$VS = $'+'{0:0.3}'.format(vs)+r', $\mu = $'+'{0:0.2}'.format(mu)+r', $\sigma_\mu = $'+'{0:0.2}'.format(sigma))

VG.set_ylabel(r'Mean $VG$')
tbar.set_ylabel(r'Mean trait value relative to optimum ($\bar{P}/P_{opt}$)')
sfs.set_xlabel('Derived mutation count')
sfs.set_ylabel('Mean proportion of mutations')
    
for ax in [VG,tbar]:
    ax.set_xlabel("Time (generations)")

VG.set_xlim(9000,13000)
tbar.set_xlim(9000,13000)
VG.legend(loc="upper left",fontsize=8)
tbar.legend(loc="lower right",fontsize=8)
sfs.legend(loc="upper right",fontsize=8)
plt.tight_layout()
plt.savefig("vary_mu_VMconstant.png")

popstats_output.close()
sfs_output.close()
