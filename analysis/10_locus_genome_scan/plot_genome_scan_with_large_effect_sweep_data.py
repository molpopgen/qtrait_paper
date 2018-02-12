import seaborn as sns
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import scipy.stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

df = pd.read_csv('raw_fixations.txt.gz',sep=' ')
df = df[df.sweep_type != 'none']
df['scaled_time'] = (df.g - 5e4)/5e3
df['ghat'] = 2.0*np.sqrt(2.0)*np.sqrt(df.mu)
#df = df[df.s.abs() > df.ghat]

choices = ['SL','SS','HL','HS']
conditions = [(df.sweep_type=='soft')&(df.s.abs()>df.ghat),
    (df.sweep_type=='soft')&(df.s.abs()<=df.ghat),
    (df.sweep_type=='hard')&(df.s.abs()>df.ghat),
    (df.sweep_type=='hard')&(df.s.abs()<=df.ghat)]
df['detailed_type']=np.select(conditions,choices)
soft=df[df.sweep_type=='soft']

df = df[(df.detailed_type=='HL')|(df.detailed_type=='SL')]

colormap = cm.viridis
#colormap = cm.PiYG
#colormap = cm.Spectral
normalize = mcolors.Normalize(vmin=0,vmax=1)
# In order to make the final colorbar meaningful,
# we need a custom palette per panel.
# We do that by mapping values for the 
# group's hue axis (s) to the normalized
# inputs to the colormap.
p = {i:colormap(normalize(scipy.stats.norm.cdf(i,scale=0.25))) for i in df.s}
    
#fig = plt.figure(figsize=(20,20))
# gs = gridspec.GridSpec(1,1)
# gs.update(left=0.1, right=0.8, wspace=0.05)
# axes = [plt.subplot(i) for i in gs]
fp = sns.factorplot(x='sweep_type',y='scaled_time',
        col='mu',row='opt',data=df,hue='s',palette=p,
        row_order=[1.0,0.5,0.1],
        #row_order=[i for i in df.opt.unique()].reverse(),
        kind='swarm',legend=False)
plt.subplots_adjust(left=0.1, right=0.8, wspace=0.05)
for ax in fp.axes.flatten():
    ax.set_ylim(-0.5,1.5)
fp.set_ylabels('Origin time of fixation')
fp.set_xlabels('Type of sweep')
#fp.savefig('test3.png')
#axes[0] = fp.axes[0][0].get_figure()
cbar_ax = fp.fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = matplotlib.colorbar.ColorbarBase(cbar_ax,cmap=colormap,
        norm=normalize)
cb.set_ticks([0.1,0.25,0.5,0.75,0.9])
cb.set_ticklabels(['{0:0.2f}'.format(scipy.stats.norm.ppf(i,scale=0.25)) for i in [0.1,0.25,0.5,0.75,0.9]])

fp.savefig('test2.png')

