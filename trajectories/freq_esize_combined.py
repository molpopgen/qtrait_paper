import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
import pandas as pd
import sqlite3
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import scipy.stats
import sys


XLIM=(np.log10(1. / (1e4)), 0.0)
YLIM=(-0.25,0.8)
mpl.rcParams.update({'font.size': 22})
fig = plt.figure(figsize=(20, 20))
files = sorted(glob.glob("*_traj_merged_recent_fixations.db"))
axes=[]
gs = gridspec.GridSpec(6, 6, width_ratios=[
                       4, 1, 4, 1, 4, 1], height_ratios=[1, 4, 1, 4, 1, 4],
                       wspace=0.0, hspace=0.0)
I = 2
J = 0

for fi, i in zip(files, range(len(files))):
    s = fi.split('_')
    opt = float(s[3])
    mu = float(s[5])
    conn = sqlite3.connect(fi)
    x = pd.read_sql(
        'select * from freqs where origin < 50000 and generation >= 50000', conn)
    g = x.groupby(['repid', 'origin', 'pos', 'esize'])[
        'generation'].transform('max').reset_index(name='maxgen')
    x['maxgen'] = g.maxgen
    x = x[x.generation == 50000]
    print(opt, mu, x.esize.min(), x.esize.max())
    x.drop_duplicates(
        subset=['repid', 'origin', 'pos', 'generation'], inplace=True)
    x.freq = x.freq.apply('log10')
    conn.close()
    print(scipy.stats.mstats.mquantiles(x.esize, prob=[1e-2, 0.5, 1. - 1e-2]))
    # Try a mplotlib kde
    ax = plt.subplot(gs[2 * I + 1, 2 * J])
    axes.append(ax)
    cs = ax.hexbin(x.freq,x.esize,gridsize=25,cmap='Blues')
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    cbaxes = inset_axes(ax, width="45%", height="3%", loc=1)
    cbar = plt.colorbar(cs, cax=cbaxes, orientation='horizontal')
    cbar.outline.set_edgecolor('black')
    cbar.outline.set_linewidth(2)
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),
                            rotation=90, color='black')
    ax.axhline(y=2 * np.sqrt(2.0) * np.sqrt(mu), linestyle='dashed',color='royalblue')
    ax.axhline(y=-2 * np.sqrt(2.0) * np.sqrt(mu), linestyle='dashed',color='royalblue')
    ax.axvline(x=np.log10(0.05), linestyle='dashed',color='royalblue')
    if mu < 1e-3:
        ax.text(-3.75, 0.6, r'$z_o = $' +
                '{0:0.2f}'.format(opt) + '\n' + r'$\mu = $' +
                '{0:0.5f}'.format(mu),
                color='black')
    else:
        ax.text(-3.75, 0.6, r'$z_o = $' +
                '{0:0.2f}'.format(opt) + '\n' + r'$\mu = $' +
                '{0:0.3f}'.format(mu),
                color='black')
    # Plot the marginals
    ax = plt.subplot(gs[2 * I, 2 * J])
    n, b, p = ax.hist(x.freq, 50, color='lightsteelblue')
    ax.set_xlim(*XLIM)
    ax.axis('off')
    for k in ax.get_xticklabels():
        k.set_visible(False)
        k.set_fontsize(0.0)
    for k in ax.get_yticklabels():
        k.set_visible(False)
        k.set_fontsize(0.0)

    ax = plt.subplot(gs[2 * I + 1, 2 * J + 1])
    n, b, p = ax.hist(x.esize, 50, orientation='horizontal',color='lightsteelblue')
    ax.set_ylim(*YLIM)
    ax.axis('off')
    for k in ax.get_xticklabels():
        k.set_visible(False)
        k.set_fontsize(0.0)
    for k in ax.get_yticklabels():
        k.set_visible(False)
        k.set_fontsize(0.0)
    J += 1
    if J == 3:
        J = 0
        I -= 1
    # p = sns.jointplot('freq', 'esize', data=x, kind='hex',
    #                   stat_func=None, xlim=(np.log10(1. / (1e4)), 0.0), ylim=(-0.25, 0.8),
    #                   space=0)
    # plt.subplots_adjust(right=0.9)
    # cbaxes = inset_axes(p.ax_joint, width="45%", height="3%", loc=1)
    # cbar = plt.colorbar(cax=cbaxes, orientation='horizontal')
    # cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=90)
    # p.set_axis_labels(
    #     r'$log_{10}$' + ' frequency at optimum shift', 'Effect size (' + r'$\gamma$' + ')')
    # p.ax_joint.text(-3.75, 0.6, r'$z_o = $' +
    #                 '{0:0.2f}'.format(opt) + '\n' + r'$\mu = $' + '{0:0.5f}'.format(mu))
    # p.ax_joint.axhline(y=2 * np.sqrt(2.0) * np.sqrt(mu), linestyle='dashed')
    # p.ax_joint.axhline(y=-2 * np.sqrt(2.0) * np.sqrt(mu), linestyle='dashed')
    # p.ax_joint.axvline(x=np.log10(0.05), linestyle='dashed')
    # # Depending on where this plot will end up when merged,
    # # we want to suppress axis label and tick output
    # if opt != 0.1:
    #     p.ax_joint.get_xaxis().set_visible(False)
    # if mu != 2.5e-4:
    #     p.ax_joint.get_yaxis().set_visible(False)
    # ofname = 'FreqVsEsizeJoint' + str(i) + '.svg'
    # p.savefig(ofname, format="svg")

for i in range(1, 9, 3):
    for j in range(i, i + 2):
        ax = axes[j]
        if j > 2:
            for k in ax.get_xticklabels():
                k.set_visible(False)
                k.set_fontsize(0.0)
        for k in ax.get_yticklabels():
            k.set_visible(False)
            k.set_fontsize(0.0)
for i in [3, 6]:
    ax = axes[i]
    for k in ax.get_xticklabels():
        k.set_visible(False)
        k.set_fontsize(0.0)

# "Finally", set axis labels for 
# external panels
    # p.set_axis_labels(
    #     r'$log_{10}$' + ' frequency at optimum shift', 'Effect size (' + r'$\gamma$' + ')')
for i in range(3):
    ax = axes[i]
    ax.set_xlabel(r'$log_{10}$' + ' frequency at optimum shift')
for i in range(0,7,3):
    ax = axes[i]
    ax.set_ylabel('Effect size (' + r'$\gamma$'+')')

plt.tight_layout()
plt.savefig('FreqVsEffectSizeStandingVariation.pdf')
