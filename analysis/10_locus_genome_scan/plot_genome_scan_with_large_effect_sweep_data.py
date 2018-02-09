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
df['scaled_time'] = (df.g - 5e4)/5e4
df['ghat'] = 2.0*np.sqrt(2.0)*np.sqrt(df.mu)
#df = df[df.s.abs() > df.ghat]

choices = ['SL','SS','HL','HS']
conditions = [(df.sweep_type=='soft')&(df.s.abs()>df.ghat),
    (df.sweep_type=='soft')&(df.s.abs()<=df.ghat),
    (df.sweep_type=='hard')&(df.s.abs()>df.ghat),
    (df.sweep_type=='hard')&(df.s.abs()<=df.ghat)]
df['detailed_type']=np.select(conditions,choices)
print(df.head())
#sys.exit(0)
soft=df[df.sweep_type=='soft']
print(soft.g.min(),soft.g.max())
# df2=df[(df.opt == 1.0)&(df.mu==2.5e-4)]
# 
# plot = sns.stripplot(x='sweep_type',y='scaled_time',hue='s',data=df2,
#         palette='viridis',jitter=True, edgecolor='none', alpha=.60)
# plot.get_legend().set_visible(False)
# sns.despine()
# normalize = mcolors.Normalize(vmin=df['s'].min(), vmax=df['s'].max())
# colormap = cm.viridis
# for n in df['s']:
#     plt.plot(color=colormap(normalize(n)))
# scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
# scalarmappaple.set_array(df['s'])
# plt.colorbar(scalarmappaple)
# fig=plot.get_figure()
# plot.set_ylabel("Time since optimum shift (units of N generations)")
# fig.savefig('test.png')

# g = sns.FacetGrid(df,col="mu",row="opt")
# g = g.map(sns.stripplot,x='sweep_type',y='scaled_time',hue='s',data=df,jitter=True)
# g.savefig('test2.png')


#DATA = df[df.opt > 0.1]

fig = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(3,3,height_ratios=(1,1,1))
gs.update(left=0.1, right=0.8, wspace=0.05)
axes = [plt.subplot(i) for i in gs]
for i in axes[1:]:
    axes[0].get_shared_y_axes().join(axes[0],i)
groups = df.groupby(['opt','mu'])
colormap = cm.viridis
normalize = mcolors.Normalize(vmin=df['s'].min(), vmax=df['s'].max())
# for n in DATA['s']:
#     plt.plot(color=colormap(normalize(n)))
# scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
# scalarmappaple.set_array(DATA['s'])
print(len(axes))
i=0
for n,g in groups:
    print(n,axes[i])
    #sns.violinplot(x='sweep_type',y='scaled_time',data=g,ax=axes[i],inner=None)
    sns.swarmplot(x='detailed_type',y='scaled_time',hue='s',data=g,
            ax=axes[i],size=3,
            #order=['soft','hard'],
            order=choices,
            palette='viridis',#jitter=True,
            edgecolor='none', alpha=1)
    axes[i].get_legend().set_visible(False)
    axes[i].set_title(r'$z_o = $'+'{0:0.2f}'.format(n[0])+', $\mu = $'+'{0:0.3f}'.format(n[1]))
    i+=1
for i in [0,1,2,3,4,5]:
    axes[i].get_xaxis().set_visible(False)
for i in [1,2,4,5,7,8]:
    axes[i].get_yaxis().set_visible(False)
axes[0].set_ylabel("Time since optimum shift")
axes[3].set_ylabel("Time since optimum shift")
axes[6].set_ylabel("Time since optimum shift")
axes[8].set_xlabel('')
axes[6].set_xlabel('')
axes[7].set_xlabel("Type of sweep")
# plt.tight_layout()
# cax,kw = matplotlib.colorbar.make_axes([ax for ax in [axes[2],axes[5]]])
# plt.colorbar(scalarmappaple,cax=cax,**kw)
#plt.colorbar(scalarmappaple,ax=[axes[2],axes[5],axes[8]])
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = matplotlib.colorbar.ColorbarBase(cbar_ax,cmap=colormap,
        norm=normalize)
# cb.set_label(r'$\gamma$',rotation='vertical')
#fig.colorbar(scalarmappaple, cax=cbar_ax)
fig.savefig('test2.png')


df
