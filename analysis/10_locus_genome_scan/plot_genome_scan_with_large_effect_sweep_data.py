import seaborn as sns
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import numpy as np
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

fig = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(3,3,height_ratios=(1,1,1))
gs.update(left=0.1, right=0.8, wspace=0.05)
axes = [plt.subplot(i) for i in gs]
for i in axes[1:]:
    axes[0].get_shared_y_axes().join(axes[0],i)
groups = df.groupby(['opt','mu'])
colormap = cm.viridis
normalize = mcolors.Normalize(vmin=df['s'].min(), vmax=df['s'].max())

print(len(axes))
i=0
for n,g in groups:
    print(n,axes[i])
    # In order to make the final colorbar meaningful,
    # we need a custom palette per panel.
    # We do that by mapping values for the 
    # group's hue axis (s) to the normalized
    # inputs to the colormap.
    p = {i:colormap(normalize(i)) for i in g.s}
    ax = sns.swarmplot(x='detailed_type',y='scaled_time',hue='s',data=g,
                ax=axes[i],size=3,
                order=choices,
                palette=p,
                edgecolor='none', alpha=1)
    axes[i].get_legend().set_visible(False)
    axes[i].set_title(r'$z_o = $'+'{0:0.2f}'.format(n[0])+', $\mu = $'+'{0:0.5f}'.format(n[1]))
    i+=1
for i in [0,1,2,3,4,5]:
    axes[i].get_xaxis().set_visible(False)
for i in [1,2,4,5,7,8]:
    axes[i].get_yaxis().set_visible(False)
axes[0].set_ylabel("Origin time of fixation")
axes[3].set_ylabel("Origin time of fixation")
axes[6].set_ylabel("Origin time of fixation")
axes[8].set_xlabel('')
axes[6].set_xlabel('')
axes[7].set_xlabel("Type of sweep")

cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = matplotlib.colorbar.ColorbarBase(cbar_ax,cmap=colormap,
        norm=normalize)

fig.savefig('test2.png')
