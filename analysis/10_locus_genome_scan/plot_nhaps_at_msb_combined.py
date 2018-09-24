import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sqlite3
import glob
import re
import sys

XLIM=(-0.25,0.8)
YLIM=(0,4100)
sns.set(style='white', font_scale=1.5)
JOIN_FIELDS = ['repid', 'opt', 'mu', 'esize', 'pos', 'origin']
fixations = pd.read_csv('raw_fixations.txt.gz', sep=' ')
fixations = fixations[(fixations.g < 50000) & (fixations.ftime > 50000)]
fixations.rename(columns={'s': 'esize', 'g': 'origin'}, inplace=True)
fixations.set_index(JOIN_FIELDS, inplace=True)
files = sorted(glob.glob('../../mlocus_pickle/*.nhaps_at_origin.sqlite3'))
print(fixations.head())
print(files)

matplotlib.rcParams.update({'font.size': 22})
fig = plt.figure(figsize=(20, 20))
gs = gridspec.GridSpec(6, 6, width_ratios=[
                       4, 1, 4, 1, 4, 1], height_ratios=[1, 4, 1, 4, 1, 4],
                       wspace=0.0, hspace=0.0)
I = 2
J = 0
dfs = []
for fi in files:
    print(fi)
    s = fi.split(".")
    opt_re = re.search('opt\d\.?\d+\.', fi)
    opt = float(fi[opt_re.start():(opt_re.end() - 1)].replace('opt', ''))
    mu_re = re.search('mu\d\.?\d+\.', fi)
    mu = float(fi[mu_re.start():(mu_re.end() - 1)].replace('mu', ''))
    c = sqlite3.connect(fi)
    df = pd.read_sql('select * from data', c)  # group by repid,pos,origin',c)
    df.drop_duplicates(inplace=True)
    df['mu'] = [mu] * len(df.index)
    df['opt'] = [opt] * len(df.index)
    dfs.append(df)
    c.close()

df = pd.concat(dfs)
df_with_fixations = df.set_index(JOIN_FIELDS)
df_with_fixations = df_with_fixations.join(fixations).dropna().reset_index()
print(df_with_fixations.nhaps.min())

gr = df_with_fixations.groupby(['mu', 'opt'])
axes = []

#for n, g in gr:
for mu in [5e-3,1e-3,2.5e-4]:
    for opt in [1.0,0.5,0.1]:
        g = gr.get_group((mu,opt))
        n = (mu,opt)
        print("processing", n, len(g.index))
        ax = plt.subplot(gs[2 * I + 1, 2 * J])
        axes.append(ax)
        cs = ax.hexbin(g.esize,g.nhaps,gridsize=25,cmap='Blues',mincnt=0)
        ax.set_xlim(*XLIM)
        ax.set_ylim(*YLIM)
        cbaxes = inset_axes(ax, width="45%", height="3%", loc=1)
        cbar = plt.colorbar(cs, cax=cbaxes, orientation='horizontal')
        cbar.outline.set_edgecolor('black')
        cbar.outline.set_linewidth(2)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),
                                rotation=90, color='black')
        if n[0] < 1e-3:
            ax.text(0.35, 2800, r'$z_o = $' +
                    '{0:0.2f}'.format(n[1]) + '\n' + r'$\mu = $' +
                    '{0:0.5f}'.format(n[0]),
                    color='black')
        else:
            ax.text(0.35, 2800, r'$z_o = $' +
                    '{0:0.2f}'.format(n[1]) + '\n' + r'$\mu = $' +
                    '{0:0.3f}'.format(n[0]),
                    color='black')
        ax.axvline(x=2 * np.sqrt(2.0) * np.sqrt(mu), linestyle='dashed',color='royalblue',
                linewidth=3)
        # Plot the marginals
        ax = plt.subplot(gs[2 * I, 2 * J])
        n, b, p = ax.hist(g.esize, 50, color='lightsteelblue')
        ax.set_xlim(*XLIM)
        ax.axis('off')
        for k in ax.get_xticklabels():
            k.set_visible(False)
            k.set_fontsize(0.0)
        for k in ax.get_yticklabels():
            k.set_visible(False)
            k.set_fontsize(0.0)

        ax = plt.subplot(gs[2 * I + 1, 2 * J + 1])
        n, b, p = ax.hist(g.nhaps, 50, orientation='horizontal',color='lightsteelblue')
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
        
    # p = sns.jointplot('esize', 'nhaps', data=g, kind='hex',
    #                   stat_func=None, xlim=(-0.25, 1), ylim=(0, 4100),
    #                   space=0)  # ,joint_kws=dict(bins=len(g.index)))
    # p.ax_joint.text(0.4, 3000, r'$\mu = $' +
    #                 '{0:0.5f}\n'.format(n[0]) + r'$z_o = $' + '{0:0.5f}'.format(n[1]))
    # p.ax_joint.set_yticks([1] + [i for i in range(500, 4000, 500)])
    # p.ax_joint.axvline(x=2 * np.sqrt(2.0) * np.sqrt(n[0]), linestyle='dashed')
    # plt.subplots_adjust(right=0.9)
    # cbaxes = inset_axes(p.ax_joint, width="45%", height="3%", loc=1)
    # cbar = plt.colorbar(cax=cbaxes, orientation='horizontal')
    # cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=90)
    # xlab = "Effect size (" + r'$\gamma$' + ")"
    # ylab = "Number of haplotypes with mutation."
    # if n[1] < 1.0:
    #     ylab = ""
    #     p.ax_joint.get_yaxis().set_visible(False)
    # if n[0] < 0.005:
    #     xlab = ""
    #     p.ax_joint.get_xaxis().set_visible(False)
    # p.set_axis_labels(xlab, ylab)

    # ofname = "NhapsVsEsizeFixationsAtMSB" + str(I) + ".svg"
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
    ax.set_xlabel('Effect size (' + r'$\gamma$'+')')
for i in range(0,7,3):
    ax = axes[i]
    ax.set_ylabel('Number of haplotypes with mutation')
plt.tight_layout()
plt.savefig("NhapsVsEsizeFixationsAtMSB.pdf")
