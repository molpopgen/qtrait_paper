import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x=pd.read_hdf('validation_summaries.h5')

sigmu = sorted(x.sigmu.unique())
sige = sorted(x.sige.unique())

linecolors = [ plt.cm.viridis(f) for f in np.linspace(0.1,0.9,4) ]
f,ax=plt.subplots(3,2,sharex='col',sharey='row')

for column in range(len(sige)):
    for row in range(len(sigmu)):
        xi=x[(x.sige==sige[column])&(x.sigmu==sigmu[row])]
        xig=xi.groupby(['r'])
        ax[row,column].plot(sorted(xi.mu.unique()),sorted(xi.EVG.unique()),ls='dotted',color='black')
        coloridx=0
        for n,g in xig:
            print n
            ax[row,column].plot(g.mu,g.meanVG,ls='solid',color=linecolors[coloridx])
            ax[row,column].set_title(r'$\sigma_E = $'+'{0:0.2f}'.format(sige[column]) +r', $\sigma_\mu = $'+'{0:0.2f}'.format(sigmu[row]))
            coloridx+=1

plt.savefig('foo.tiff')
