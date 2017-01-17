import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x=pd.read_hdf('validation_summaries.h5')

sigmu = sorted(x.sigmu.unique())

linecolors = [ plt.cm.viridis(f) for f in np.linspace(0.1,0.9,4) ]
f,ax=plt.subplots(3,sharex='col',sharey='row')

for row in range(len(sigmu)):
    xi=x[(x.sigmu==sigmu[row])]
    xig=xi.groupby(['r'])
    ax[row].plot(sorted(xi.mu.unique()),sorted(xi.EVG.unique()),ls='dotted',color='black')
    coloridx=0
    for n,g in xig:
        print n
        ax[row].plot(g.mu,g.meanVG,ls='solid',color=linecolors[coloridx],label=r'$r = $'+'{0:0.2f}'.format(n))
        ax[row].set_title(r'$\sigma_\mu = $'+'{0:0.3f}'.format(sigmu[row]))
        coloridx+=1

for axi in range(len(sigmu)):
    ax[axi].legend(loc='best',prop={'size':10})

plt.savefig('foo.tiff')
