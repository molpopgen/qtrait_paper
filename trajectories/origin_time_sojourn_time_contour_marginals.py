import matplotlib as mpl
mpl.use('Agg')
import pandas as pd
import sqlite3
import glob
import matplotlib.mlab
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import scipy.stats
import sys


mpl.rcParams.update({'font.size': 22})

files = sorted(glob.glob("*_traj_merged_recent_fixations.db"))
fig = plt.figure(figsize=(20, 20))
# gs = gridspec.GridSpec(3,6,width_ratios=[4,1,4,1,4,1])
gs = gridspec.GridSpec(6, 6, width_ratios=[
                       4, 1, 4, 1, 4, 1], height_ratios=[1, 4, 1, 4, 1, 4],
                       wspace=0.0, hspace=0.0)
I = 2
J = 0

# store the "main" axes:
axes = []
for fi, i in zip(files, range(len(files))):
    s = fi.split('_')
    opt = float(s[3])
    mu = float(s[5])
    conn = sqlite3.connect(fi)
    x = pd.read_sql(
        'select * from freqs where freq == 1.0', conn)
    x['sojourn_time'] = (x.generation - x.origin + 1.0) / 5e3
    x['scaled_origin'] = (x.origin - 5e4)
    print(opt, mu, x.esize.min(), x.esize.max())
    conn.close()

    # Try a mplotlib kde
    ax = plt.subplot(gs[2 * I + 1, 2 * J])
    axes.append(ax)
    cs = ax.hexbin(x.scaled_origin,x.sojourn_time,gridsize=25,cmap='Blues')
    ax.set_xlim(-4*5e3,5e3*4.01)
    ax.set_ylim(0,8.01)
    cbaxes = inset_axes(ax, width="45%", height="3%", loc=1)
    cbar = plt.colorbar(cs, cax=cbaxes, orientation='horizontal')
    cbar.outline.set_edgecolor('black')
    cbar.outline.set_linewidth(2)
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),
                            rotation=45, color='black')
    if mu < 1e-3:
        ax.text(-3.5, 6.5, r'$z_o = $' +
                '{0:0.2f}'.format(opt) + '\n' + r'$\mu = $' +
                '{0:0.5f}'.format(mu),
                color='black')
    else:
        ax.text(-3.5, 6.5, r'$z_o = $' +
                '{0:0.2f}'.format(opt) + '\n' + r'$\mu = $' +
                '{0:0.3f}'.format(mu),
                color='black')
    # Plot the marginals
    ax = plt.subplot(gs[2 * I, 2 * J])
    n, b, p = ax.hist(x.scaled_origin, 50, color='lightsteelblue')
    ax.set_xlim(-4*5e3, 5e3*4)
    ax.axis('off')
    for k in ax.get_xticklabels():
        k.set_visible(False)
        k.set_fontsize(0.0)
    for k in ax.get_yticklabels():
        k.set_visible(False)
        k.set_fontsize(0.0)

    ax = plt.subplot(gs[2 * I + 1, 2 * J + 1])
    n, b, p = ax.hist(x.sojourn_time, 50, orientation='horizontal',color='lightsteelblue')
    ax.set_ylim(0, 8)
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

# Hid axes on "internal" panels
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
for i in range(3):
    ax = axes[i]
    for t in ax.get_xticklabels():
        t.set_rotation(45)
axes[1].set_xlabel('Origin time of fixation\n(generations since optimum shift)')
for i in range(0,7,3):
    ax = axes[i]
axes[3].set_ylabel('Fixation time, in units of N generations')


plt.tight_layout()
plt.savefig('OriginTimeSojournTime.pdf')
