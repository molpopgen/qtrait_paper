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


files = glob.glob("*_traj_merged_recent_fixations.db")

for fi, i in zip(files, range(len(files))):
    s = fi.split('_')
    opt = float(s[3])
    mu = float(s[5])
    conn = sqlite3.connect(fi)
    x = pd.read_sql(
        'select * from freqs where freq == 1.0', conn)
    x['sojourn_time'] = (x.generation-x.origin + 1.0)/5e3
    x['scaled_origin'] = (x.origin - 5e4)/5e3
    print(opt, mu, x.esize.min(), x.esize.max())
    conn.close()
    fig = plt.figure(figsize=(10, 10))
    p = sns.jointplot('scaled_origin','sojourn_time', data=x, kind='hex',
                      stat_func=None, xlim=(-4,4),ylim=(0,8),
                      space=0)
    plt.subplots_adjust(right=0.9)
    cbaxes = inset_axes(p.ax_joint, width="45%", height="3%", loc=1)
    cbar = plt.colorbar(cax=cbaxes, orientation='horizontal')
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=90)
    p.set_axis_labels("Origin time of fixation (units of N generations)","Time to fixation (units of N generations)")
    # p.set_axis_labels(
    #     r'$log_{10}$' + ' frequency at optimum shift', 'Effect size (' + r'$\gamma$' + ')')
    p.ax_joint.text(-3.5, 7, r'$z_o = $' +
                    '{0:0.2f}'.format(opt) + '\n' + r'$\mu = $' + '{0:0.5f}'.format(mu))
    p.ax_marg_y.axhline(y=x.sojourn_time.mean(), linestyle='dashed',linewidth=4)
    p.ax_marg_y.axhline(y=x.sojourn_time.median(), linestyle='dotted',linewidth=4)
    # p.ax_joint.axvline(x=np.log10(0.05), linestyle='dashed')
    # Depending on where this plot will end up when merged,
    # we want to suppress axis label and tick output
    if opt != 0.1:
        p.ax_joint.get_xaxis().set_visible(False)
    if mu != 2.5e-4:
        p.ax_joint.get_yaxis().set_visible(False)
    ofname = 'OriginVsSojournJoint' + str(i) + '.svg'
    p.savefig(ofname, format="svg")
    #sys.exit(0)

