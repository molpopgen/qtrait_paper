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

for statfile, trajfile, ofile in zip(popstatFiles, trajFiles, oFiles):
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
    matplotlib.rcParams.update({'font.size': 18})
    matplotlib.rcParams.update({'legend.fontsize': 14})
    ##Create a colormap
    #cmap = matplotlib.cm.get_cmap('Oranges')
    cmap = matplotlib.cm.get_cmap('viridis')


    ##Create a grid: top row will be pop-level parameters, bottom mutation frequencies
    from matplotlib import gridspec
    fig=plt.figure(figsize=(20,10))
    gs = gridspec.GridSpec(3,1,height_ratios=(1,2,2))

    axTOP = plt.subplot(gs[0])
    #PLot mean trait value
    axTOP.plot(pstats.scaled_time[pstats.stat=='tbar'],
             pstats.value[pstats.stat=='tbar'],
                 color='blue',alpha=0.75,linewidth=2.5,label="Mean trait value")
    #axTOP.plot(pstats.scaled_time[pstats.stat=='wbar'],
    #         pstats.value[pstats.stat=='wbar'],
    #             color='green',alpha=0.75,linewidth=2.5,label="Mean fitness")
    axTOP.plot(pstats.scaled_time[pstats.stat=='varw'],
             pstats.value[pstats.stat=='varw'].multiply(100.0),
                 color='green',alpha=0.75,linewidth=2.5,label=r'$100 \times V(w)$')
    #10xVG, so that it shows up
    axTOP.plot(pstats.scaled_time[pstats.stat=='tbar'],
             10.0*pstats.value[pstats.stat=='VG'],
                 color='purple',linestyle='solid',linewidth=2.5,alpha=1,label=r'$10 \times V(G)$')
    plt.ylabel("Value")
    plt.legend(frameon=False)#loc='upper left')
    axBOTTOM = plt.subplot(gs[1],sharex=axTOP)

    ##SORT FIXATIONS BY EFFECT SIZE
    FIXATIONS=[]
    MAX_ESIZE = 0
    for i in fixations.groups:
        fix_i = fixations.get_group(i)
        MAX_ESIZE=max(MAX_ESIZE,math.fabs(fix_i.esize.max()))
        FIXATIONS.append(fix_i)
    FIXATIONS=sorted(FIXATIONS,key = lambda F : F.esize.min(),reverse=True)
    #plot fixations
    FSTYLES=['solid','dashed','solid']  #hand-coded hack

    ENTRIES = 0
    for fix_i in FIXATIONS:
        FSTYLE=0
        #fix_i=fixations.get_group(i)
        fix_color="red"
        fix_alpha=0.25
        fix_style=FSTYLES[FSTYLE]
        if fix_i.origin.min()>50000:
            fix_color='blue'
            fix_alpha=1.0
            fix_style='solid'
        else:
            FSTYLE+=1
        esize=fix_i.esize.mean()
        origin=fix_i.origin.mean()
        ##Do a super-dirty calculation of s:
        ftime = len(fix_i.freq)
        s = -1.0*math.log(1.0/10000.0)/(float(ftime))
        #print(s)

        label = "_nolabel_"
        if math.fabs(esize) > 0.1:
             label=r'$\gamma = $'+'{0:.2f}'.format(esize)+ r', $o = $'+'{0:0.4f}'.format((origin-5e4)/5e3)
        axBOTTOM.plot(fix_i.scaled_time,fix_i.freq,#color=fix_color,
                     alpha=min(1.0,4.0*math.fabs(esize)),
                     #linestyle=fix_style,
                     linewidth=2,
                     color = cmap(1.0 - math.fabs(esize)/MAX_ESIZE),
                     #color=cmap(100.0*(esize*esize)),
                     label=label)
        # "Fixed: "+r'$e = $'+'{0:.2f}'.format(esize)+
        #              r', $o = $'+'{0:0.2f}'.format((origin-5e4)/5e3))
        #     axBOTTOM.plot(fix_i.scaled_time,fix_i.freq,#color=fix_color,
        #                  alpha=min(1.0,4.0*math.fabs(esize)),
        #                  #linestyle=fix_style,
        #                  linewidth=2,
        #                  color = cmap(1.0 - math.fabs(esize)/MAX_ESIZE),
        #                  label = "_nolegend_")
    axBOTTOM.legend(loc="upper right",ncol=1,title=r'$\gamma = \mathrm{effect\ size}$,'+'\n'+ r'$o = \mathrm{origination\ time}$',fontsize='small',frameon=False)
    #axBOTTOM.legend(bbox_to_anchor=(1.2,1),ncol=1,title=r'$e = \mathrm{effect\ size}$,'+'\n'+ r'$o = \mathrm{origination\ time}$',fontsize='small')
    axBOTTOM.set_ylabel("Mutation frequency")
        
    #Now, deal with losses
    #We will collect them into a list, for sorting
    #I really wish I could just sort the groups...
    axBOTTOM = plt.subplot(gs[2],sharex=axTOP)
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
        if math.fabs(esize) > 0.1 and ENTRIES < 5:
            label=r'$\gamma = $'+'{0:0.2f}'.format(esize)+r', $o = $'+'{0:0.4f}'.format((origin-5e4)/5e3)
            ENTRIES += 1
        axBOTTOM.plot(i.scaled_time,i.freq,
                     linewidth=2,
                     color = cmap(1.0 - math.fabs(esize)/MAX_ESIZE),
                     alpha=min(1.0,4.0*math.fabs(esize)),label=label)
        
    plt.xlim(-0.1,1)

    plt.ylim(-0.05,1)
    # axBOTTOM.text(-0.4,0.9,'Soft sweep,\n'+r'$e = 0.12,$'+'\n'+r'$o = -0.51$')#+',\n'+r'$s\approx 0.01$')
    # axBOTTOM.text(0.375,0.9,'Hard sweep,\n'+r'$e = 0.09,$'+'\n'+r'$o = 0.07$')#+',\n'+r'$s\approx 0.008$')
    # axBOTTOM.text(3.3,0.9,'Soft sweep,\n'+r'$e = 0.06$,'+'\n'+r'$o = -0.23$')#+',\n'+r'$s\approx 0.002$')
    plt.xlabel("Time since optimum shift (units of N generations)")
    axBOTTOM.set_ylabel("Mutation frequency")
    axBOTTOM.legend(loc="upper right",ncol=1,title=r'$\gamma = \mathrm{effect\ size}$,'+'\n'+ r'$o = \mathrm{origination\ time}$',fontsize='small',frameon=False)
    # axBOTTOM.hlines(0.35,0.75,3.5)
    # axBOTTOM.vlines(0.75,0.325,0.35)
    # axBOTTOM.vlines(3.5,0.325,0.35)
    # axBOTTOM.text(0.75,0.365,"After optimum is reached, "+
    #              "there is a period of long-lived mutations\naffecting fitness. "+
    #              "Genetic load is getting sorted out.",fontsize=14)
    # axTOP.vlines(3.3,0.475,0.525)
    # #ax=plt.axes()
    # axTOP.arrow(3.3,0.5,1.5,0)
    # axTOP.text(3.325,0.2,"Optimum slightly over-shot,\n"+r'$G=2\sum_i e_i = 0.54$'+' due to sweeps.')

    plt.tight_layout()

    ##Adjusts spacing b/w axes
    fig.subplots_adjust(hspace=0)   
    for ax in [axTOP]:
        plt.setp(ax.get_xticklabels(), visible=False)
        # The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
        ax.set_yticks(ax.get_yticks()[1:])  
        
    #axBOTTOM.legend(bbox_to_anchor=(1.2,1),ncol=1,title=r'$e = \mathrm{effect\ size}$,'+'\n'+ r'$o = \mathrm{origination\ time}$',fontsize='small')
    #plt.show()
    plt.savefig(ofile,bbox_inches='tight')

