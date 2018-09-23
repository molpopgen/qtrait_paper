import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sqlite3
import glob

sns.set(font_scale=1.5, style="white")
files = glob.glob('*_traj_merged.db')

dfs = []
for fi in files:
    s = fi.split("_")
    opt = s[3]
    if opt == '1':
        mu = float(s[5])
        print(fi,mu)
        c = sqlite3.connect(fi)
        df = pd.read_sql('select repid,freq,origin,esize from freqs where repid < 1 and generation == 50000',c) # group by repid,pos,origin',c)
        df.drop_duplicates(inplace=True)
        df['mu']=[mu]*len(df.index)
        dfs.append(df)
        c.close()

df = pd.concat(dfs)

gr = df.groupby(['mu'])

I=0
for n,g in gr:
    print(n,len(g.index))
    fig = plt.figure(figsize=(10, 10))
    p = sns.jointplot('esize','freq', data=g, kind='hex',
                      stat_func=None, xlim=(-1,1),ylim=(0,0.2),
                      space=0)#,joint_kws=dict(bins=len(g.index)))
    p.ax_joint.text(-0.9,0.185,r'$\mu = $'+'{0:0.5f}'.format(n))
    plt.subplots_adjust(right=0.9)
    cbaxes = inset_axes(p.ax_joint, width="45%", height="3%", loc=1)
    cbar = plt.colorbar(cax=cbaxes, orientation='horizontal')
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=90)
    p.set_axis_labels("Effect size ("+r'$\gamma$'+")","Mutant allele frequency")
    ofname="FreqEsizePriorToShift"+str(I)+".svg"
    p.savefig(ofname, format="svg")
    I+=1

    # fig = plt.figure(figsize=(10, 10))
    # plt.hist2d(x=g.esize,y=g.freq,normed=True,bins=50)
    # ofname="FreqEsizePriorToShift"+str(I)+".hist2d.svg"
    # plt.savefig(ofname)
# g = sns.FacetGrid(df,col='mu',margin_titles=True)
# g.map(plt.hist2d,x='esize',y='freq',bins=50,normed=True)
# ofname="FreqEsizePriorToShift"+str(I)+".facet.svg"
# p.savefig(ofname, format="svg")
