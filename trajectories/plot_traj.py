import matplotlib
matplotlib.use('Agg')
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sqlite3
import sys
import math
# Read in VG, etc.

popstatFiles = ['../sqlite/H2_1.0_OPT_1_mu_0.005_sigmu0.25_stats.db','../sqlite/H2_1.0_OPT_1_mu_0.001_sigmu0.25_stats.db','../sqlite/H2_1.0_OPT_1_mu_0.00025_sigmu0.25_stats.db']
trajFiles = ['H2_1.0_OPT_1_mu_0.005_sigmu0.25_traj_merged.db','H2_1.0_OPT_1_mu_0.001_sigmu0.25_traj_merged.db','H2_1.0_OPT_1_mu_0.00025_sigmu0.25_traj_merged.db']
oFiles = ['trajectories_himu.pdf', 'trajectories_medmu.pdf', 'trajectories_lomu.pdf']
# popstatFiles = ['../sqlite/H2_1.0_OPT_1_mu_0.00025_sigmu0.25_stats.db']
# trajFiles = ['H2_1.0_OPT_1_mu_0.00025_sigmu0.25_traj_merged.db']
# oFiles = ['H2_1.0_OPT_1_mu_0.00025_sigmu0.25_traj_trajectories.pdf']

matplotlib.rcParams.update({'font.size': 6})
matplotlib.rcParams.update({'legend.fontsize': 4})
##Create a colormap
cmap = matplotlib.cm.get_cmap('viridis')

from matplotlib import gridspec
#fig=plt.figure(figsize=(20,30))
fig=plt.figure(figsize=(12,6))
fig = plt.figure()
gs = gridspec.GridSpec(3,3,height_ratios=(0.5,1,1))
#Set up the axes here
_T1 = plt.subplot(gs[0])
_T2 = plt.subplot(gs[1],sharex=_T1,sharey=_T1)
_T3 = plt.subplot(gs[2],sharex=_T1,sharey=_T1)
_M1 = plt.subplot(gs[3],sharex=_T1)
_M2 = plt.subplot(gs[4],sharex=_T1,sharey=_M1)
_M3 = plt.subplot(gs[5],sharex=_T1,sharey=_M1)
_B1 = plt.subplot(gs[6],sharex=_T1)
_B2 = plt.subplot(gs[7],sharex=_T1,sharey=_B1)
_B3 = plt.subplot(gs[8],sharex=_T1,sharey=_B1)

# Arrays of axes for "easy" processing below
TOPS=[_T1,_T2,_T3]
MIDS=[_M1,_M2,_M3]
BOTTOMS=[_B1,_B2,_B3]

# Set xlims and ylims
_T1.set_xlim(-0.1,1)
_T1.set_ylim(-0.05,1.25)
_M1.set_ylim(-0.05,1)
_B1.set_ylim(-0.05,1)

#Set labels
_T1.set_ylabel("Value")
for ax in [_M1,_B1]:
    ax.set_ylabel("Mutation frequency")
_B2.set_xlabel("Time since optimum shift (units of N generations)")

# Disable axis ticks on internal panels
for ax in TOPS + MIDS:
    ax.get_xaxis().tick_bottom()
    ax.get_xaxis().set_visible(False)
for ax in TOPS[1:] + MIDS[1:] + BOTTOMS[1:]:
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().set_visible(False)

# Row titles
_T3.set_title(r'$\mu = 0.005, z_0 = 1, \hat\gamma = $'+'{0:0.3f}'.format(2.0*math.sqrt(2)*math.sqrt(0.005)))
_T2.set_title(r'$\mu = 0.001, z_0 = 1, \hat\gamma = $'+'{0:0.3f}'.format(2.0*math.sqrt(2)*math.sqrt(0.001)))
_T1.set_title(r'$\mu = 2.5\times 10^{-4}, z_0 = 1, \hat\gamma = $'+'{0:0.3f}'.format(2.0*math.sqrt(2)*math.sqrt(0.00025)))

# Row counters
TOP=0
MID=0
BOTTOM=0

muvals = [0.005, 0.001, 0.00025]
for statfile, trajfile,mu in zip(reversed(popstatFiles), reversed(trajFiles), reversed(muvals)):

    axTOP = TOPS[TOP]
    axMID = MIDS[MID]
    axBOTTOM = BOTTOMS[BOTTOM]

    ghat = 2.0*math.sqrt(2.0)*math.sqrt(mu)
    print("starting",ghat)

    # Read in trait value data for this rep
    con = sqlite3.connect(statfile)
    pstats=pd.read_sql('select * from data where rep == 42 and generation >= 45000 and generation <= 55000',con)
    con.close()
    pstats['scaled_time'] = pd.Series((pstats.generation - 50000),dtype=np.float64)/5000.0

    # Read in trajectories for this rep
    con = sqlite3.connect(trajfile)
    data=pd.read_sql('select * from freqs where repid == 42 and origin >= 45000 and origin <= 55000',con)
    con.close()
    data['scaled_time'] = pd.Series((data.generation - 50000),dtype=np.float64)/5000.0

    g = data.groupby(['esize','origin','pos'])

    #Filter groups into new pd.DF based on fixed vs. lost, then re-group
    fixations=g.filter(lambda x:x['freq'].max()==1.0).groupby(['esize','origin','pos'])

    #Filter any mutation that didn't live long enough or get to at least a certaing frequency, and did not fix
    losses=g.filter(lambda x:len(x['freq'])>=5 and x['freq'].max()>=0.01 and x['freq'].max()<1.0).groupby(['esize','origin','pos'])

    #PLot mean trait value
    axTOP.plot(pstats.scaled_time[pstats.stat=='tbar'],
             pstats.value[pstats.stat=='tbar'],
                 color='blue',alpha=0.75,linewidth=1.5,label="Mean trait value")
    # axTOP.plot(pstats.scaled_time[pstats.stat=='varw'],
    #          pstats.value[pstats.stat=='varw'].multiply(100.0),
    #              color='green',alpha=0.75,linewidth=2.5,label=r'$100 \times V(w)$')
    #10xVG, so that it shows up
    axTOP.plot(pstats.scaled_time[pstats.stat=='tbar'],
             5.0*pstats.value[pstats.stat=='VG'],
                 color='purple',linestyle='solid',linewidth=1.5,alpha=1,label=r'$5 \times V(G)$')
    axTOP.legend(loc='center right',frameon=False)#loc='upper left')

    ##SORT FIXATIONS BY EFFECT SIZE
    FIXATIONS=[]
    MAX_ESIZE = 0
    for i in fixations.groups:
        fix_i = fixations.get_group(i)
        MAX_ESIZE=max(MAX_ESIZE,math.fabs(fix_i.esize.max()))
        FIXATIONS.append(fix_i)
    FIXATIONS=sorted(FIXATIONS,key = lambda F : F.esize.min(),reverse=True)
    #plot fixations

    ENTRIES=0
    for fix_i in FIXATIONS:
        esize=fix_i.esize.mean()
        origin=fix_i.origin.mean()

        label = "_nolabel_"
        if math.fabs(esize) > ghat and ENTRIES < 5:
            label=r'$\gamma = $'+'{0:.2f}'.format(esize)+ r', $o = $'+'{0:0.4f}'.format((origin-5e4)/5e3)
            ENTRIES += 1
        axMID.plot(fix_i.scaled_time,fix_i.freq,#color=fix_color,
                     alpha=min(1.0,4.0*math.fabs(esize)),
                     #linestyle=fix_style,
                     linewidth=1.5,
                     color = cmap(1.0 - math.fabs(esize)/MAX_ESIZE),
                     #color=cmap(100.0*(esize*esize)),
                     label=label)
    axMID.legend(loc="center right",ncol=1,title=r'$\gamma = \mathrm{effect\ size}$,'+'\n'+ r'$o = \mathrm{origination\ time}$',fontsize='small',frameon=False)

    #Now, deal with losses
    #We will collect them into a list, for sorting
    #I really wish I could just sort the groups...
    LOSSES=[]
    for i in losses.groups:
        losses_i=losses.get_group(i)
        if losses_i.origin.max() >= 45000 and losses_i.origin.max() < 55000:
            LOSSES.append(losses_i)
            
    LOSSES=sorted(LOSSES,key = lambda LI: LI.esize.min()*LI.esize.min(),reverse=True) #Sort base on origin time, .min is again a hack
    ENTRIES=0
    for i in LOSSES:
        esize=i.esize.mean()
        origin=i.origin.mean()
        label = "_nolegend_"
        if math.fabs(esize) > ghat and ENTRIES < 5:
            label=r'$\gamma = $'+'{0:0.2f}'.format(esize)+r', $o = $'+'{0:0.4f}'.format((origin-5e4)/5e3)
            ENTRIES += 1
        axBOTTOM.plot(i.scaled_time,i.freq,
                     linewidth=1.5,
                     color = cmap(1.0 - math.fabs(esize)/MAX_ESIZE),
                     alpha=min(1.0,4.0*math.fabs(esize)),label=label)
         
    axBOTTOM.legend(loc="upper right",ncol=1,title=r'$\gamma = \mathrm{effect\ size}$,'+'\n'+ r'$o = \mathrm{origination\ time}$',fontsize='small',frameon=False)

    TOP+=1
    MID+=1
    BOTTOM+=1

gs.update(wspace=0.00, hspace=0.025)
gs.tight_layout(fig)
plt.savefig("trajectories.pdf",bbox_inches='tight')

