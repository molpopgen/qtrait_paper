import pandas as pd
import glob
import sqlite3
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.stats
import sys

files = glob.glob("*_traj_merged_recent_fixations.db")


fdata = []
for i in files:
    s = i.split('_')
    opt=float(s[3])
    mu=float(s[5])
    print(opt,mu)
    conn = sqlite3.connect(i)
    x = pd.read_sql('select *,count(repid) as sojourn_time from freqs where origin < 60000 group by repid, origin, pos',conn)
    #x.groupby(['repid','origin','pos','esize']).size().reset_index(name='sojourn_time')
    x['mu']=[mu]*len(x.index)
    x['opt']=[opt]*len(x.index)
    fdata.append(x)
    conn.close()

fixations = pd.concat(fdata)
fixations['scaled_time'] = (fixations.origin - 5e4)/5e3
fixations['scaled_sojourn_time'] = (fixations.sojourn_time)/5e3
fixations['ghat'] = 2.0*np.sqrt(2.0)*np.sqrt(fixations.mu)

choices = ['hard','soft']
conditions = [(fixations.scaled_time > 0),
        (fixations.scaled_time < 0) & (fixations.scaled_time + fixations.scaled_sojourn_time > 0.0)]
fixations['sweep_type'] = np.select(conditions,choices,default='none')
fixations = fixations[fixations.sweep_type != 'none']

choices = ['SL','SS','HL','HS']
conditions = [(fixations.sweep_type=='soft')&(fixations.esize.abs()>fixations.ghat),
    (fixations.sweep_type=='soft')&(fixations.esize.abs()<=fixations.ghat),
    (fixations.sweep_type=='hard')&(fixations.esize.abs()>fixations.ghat),
    (fixations.sweep_type=='hard')&(fixations.esize.abs()<=fixations.ghat)]
fixations['detailed_type']=np.select(conditions,choices)

# Create a color palette mapping effect sizes
# to the cdf of the DES
colormap = cm.viridis
normalize = mcolors.Normalize(vmin=0,vmax=1)
p = {i:colormap(normalize(scipy.stats.norm.cdf(i,scale=0.25))) for i in fixations.esize}
    
fixations=fixations[(fixations.detailed_type=='HL')|(fixations.detailed_type=='SL')]
fig = plt.figure(figsize=(10,10))
fp = sns.factorplot(x='sweep_type',y='scaled_time',
        col='mu',row='opt',data=fixations,hue='esize',palette=p,
        row_order=[1.0,0.5,0.1],
        #row_order=[i for i in df.opt.unique()].reverse(),
        kind='swarm',legend=False)
plt.subplots_adjust(left=0.1, right=0.8, wspace=0.05)
zo=['0.1','0.5','1.0']
mu=['0.00025','0.001','0.005']
for ax in fp.axes.flatten():
    # Here, we reformat the
    # subplot titles to use
    # the appropriate notation
    t = ax.get_title()
    zz=None
    for z in zo:
        if z in t:
            zz=z
    mm=None
    for m in mu:
        if m in t:
            mm=m
    ax.set_title(r'$z_o = $'+zz+', '+r'$\mu = $'+mm)
    ax.set_ylim(-0.5,1.5)
fp.set_ylabels('Origin time of fixation')
fp.set_xlabels('Type of sweep')
cbar_ax = fp.fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = matplotlib.colorbar.ColorbarBase(cbar_ax,cmap=colormap,
        norm=normalize)
# Set some quantiles of the DES and label
# colorbar accordingly
cb.set_ticks([0.1,0.25,0.5,0.75,0.9])
cb.set_ticklabels(['{0:0.2f}'.format(scipy.stats.norm.ppf(i,scale=0.25)) for i in [0.1,0.25,0.5,0.75,0.9]])

fp.savefig('TrajectoriesLargeEffectSweepsOriginTimes.pdf')

