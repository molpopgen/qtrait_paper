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

JOIN_FIELDS=['repid','opt','mu','esize','pos','origin']
fixations = pd.read_csv('raw_fixations.txt.gz',sep=' ')
fixations=fixations[(fixations.g < 50000) & (fixations.ftime > 50000)]
fixations.rename(columns={'s':'esize','g':'origin'},inplace=True)
fixations.set_index(JOIN_FIELDS,inplace=True)
files = glob.glob('../../mlocus_pickle/*.nhaps_at_origin.sqlite3')

dfs = []
for fi in files:
    s = fi.split(".")
    opt_re = re.search('opt\d\.?\d+\.',fi)
    opt=float(fi[opt_re.start():(opt_re.end()-1)].replace('opt',''))
    mu_re = re.search('mu\d\.?\d+\.',fi)
    mu=float(fi[mu_re.start():(mu_re.end()-1)].replace('mu',''))
    c = sqlite3.connect(fi)
    df = pd.read_sql('select * from data',c) # group by repid,pos,origin',c)
    df.drop_duplicates(inplace=True)
    df['mu']=[mu]*len(df.index)
    df['opt']=[opt]*len(df.index)
    dfs.append(df)
    c.close()

df = pd.concat(dfs)
df_with_fixations=df.set_index(JOIN_FIELDS)
df_with_fixations = df_with_fixations.join(fixations).dropna().reset_index()
print(df_with_fixations.nhaps.min())

# def add_class(row):
#     if np.isnan(row.ftime):
#         return "Loss"
#     return "Fixation"
# 
# col = df_with_fixations.apply(lambda x : add_class(x), axis = 1)
# df_with_fixations['Class'] = col
# 
# print(df_with_fixations.head())
# print(df_with_fixations[df_with_fixations['Class']=="Fixation"].head())
# 
#  
# gr = df_with_fixations.groupby(['mu','opt'])
# I=0
# gs = gridspec.GridSpec(3,3)
# column = 0
# row = 0
# plt.figure(figsize=(10,10)) 
# for n,g in gr:
#     print(n)
#     ax = plt.subplot(gs[row,column])
#     n = g.nhaps[g.Class == "Loss"]
#     e = g.esize[g.Class == "Loss"]
#     print(n.mean(),n.median())
#     ax.plot(e,n,'bo',alpha=0.1)
#     n = g.nhaps[g.Class == "Fixation"]
#     e = g.esize[g.Class == "Fixation"]
#     print(n.mean(),n.median())
#     ax.plot(e,n,'rs',alpha=1)
# 
#     #ax.hist([g.nhaps[g.Class == "Loss"],g.nhaps[g.Class=="Fixation"]],label=["Loss","Fixation"],bins=100,normed=True)
# 
#     ax.legend()
#     column += 1
#     if column % 3 == 0.0:
#         row += 1
#         column=0
# plt.savefig("foo.pdf")
# plt.close()
# I=0
# for n,g in gr:
#     print(n,len(g.index))
#     fig = plt.figure(figsize=(10, 10))
#     p = sns.jointplot('esize','nhaps', data=g, kind='hex',
#                       stat_func=None, xlim=(-1,1),
#                       space=0)#,joint_kws=dict(bins=len(g.index)))
#     #p.ax_joint.text(-0.9,0.185,r'$\mu = $'+'{0:0.5f}'.format(n))
#     plt.subplots_adjust(right=0.9)
#     cbaxes = inset_axes(p.ax_joint, width="45%", height="3%", loc=1)
#     cbar = plt.colorbar(cax=cbaxes, orientation='horizontal')
#     cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=90)
#     p.set_axis_labels("Effect size ("+r'$\gamma$'+")","Number of haplotypes with mutation.")
#     ofname="NhapsVsEsizeAtMSB"+str(I)+".svg"
#     p.savefig(ofname, format="svg")
#     I+=1
# 
gr = df_with_fixations.groupby(['mu','opt'])

I=0
for n,g in gr:
    print(n,len(g.index))
    fig = plt.figure(figsize=(10, 10))
    p = sns.jointplot('esize','nhaps', data=g, kind='hex',
                      stat_func=None, xlim=(-1,1),ylim=(0,4000),
                      space=0)#,joint_kws=dict(bins=len(g.index)))
    p.ax_joint.text(-0.9,3500,r'$\mu = $'+'{0:0.5f}\n'.format(n[0]) + r'$z_o = $' +'{0:0.5f}'.format(n[1]))
    p.ax_joint.set_yticks([1] + [i for i in range(500,4000,500)])
    plt.subplots_adjust(right=0.9)
    cbaxes = inset_axes(p.ax_joint, width="45%", height="3%", loc=1)
    cbar = plt.colorbar(cax=cbaxes, orientation='horizontal')
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=90)
    p.set_axis_labels("Effect size ("+r'$\gamma$'+")","Number of haplotypes with mutation.")
    ofname="NhapsVsEsizeFixationsAtMSB"+str(I)+".svg"
    p.savefig(ofname, format="svg")
    I+=1

