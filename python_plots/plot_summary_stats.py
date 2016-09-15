import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob,sys

x=pd.read_hdf(sys.argv[1]) #HDF5 file with mean values of summary stat per generation
STUB=sys.argv[2] #The image for a summary stat will be in file STAT+STUB+.tif
x['scaled_time']=x.generation.subtract(50000.).divide(2e4)
xg=x.groupby(['mu','opt'])

linestyles=['solid','dashed','dotted']
linecolors = [ plt.cm.viridis(f) for f in np.linspace(0.1,0.9, 3) ]
styles = dict()
colors = dict()
I=0

for i in x.mu.unique():
    colors[i]=linecolors[I]
    I+=1
I=0
for i in x.opt.unique():
    styles[i]=linestyles[I]
    I+=1

for stat in x.variable.unique():
    fig=plt.figure()
    lo=plt.subplot2grid((1,3),(0,0))
    mid=plt.subplot2grid((1,3),(0,1),sharex=lo,sharey=lo)
    hi=plt.subplot2grid((1,3),(0,2),sharex=lo,sharey=lo)
    splots=[lo,mid,hi]
    maxStat=x[(x.variable==stat) & (x.scaled_time > -0.1) & (x.scaled_time <1)].max().value
    minStat=x[(x.variable==stat) & (x.scaled_time > -0.1) & (x.scaled_time <1)].min().value
    for n,g in xg:
        SPLOT=0
        if n[1]=="0.5":
            SPLOT=1
        elif n[1]=="1":
            SPLOT=2
        splots[SPLOT].plot(g[g.variable==stat].scaled_time,
                           g[g.variable==stat].value,
                           label=r'$\mu = $'+'{0:0.1e}'.format(float(n[0])),
                           color=colors[n[0]],
                           linewidth=2)
        splots[SPLOT].legend(loc='best',fontsize=6,title=r'$P_o = $'+n[1])
    for i in splots:
        labels=i.get_xticklabels()
        plt.setp(labels, rotation=70)
    
    lo.set_xlim(-0.25,1.0)
    lo.set_ylim(minStat,maxStat)
    plt.setp(mid.get_yticklabels(),visible=False)
    plt.setp(hi.get_yticklabels(),visible=False)
    mid.set_xlabel("Time since optimum shift (units of N generations)")
    lo.set_ylabel("Mean value (per locus)")
    plt.gcf().subplots_adjust(bottom=0.225)
    fname='images/'+stat+STUB+'.tif'
    plt.savefig(fname,dpi=600)

