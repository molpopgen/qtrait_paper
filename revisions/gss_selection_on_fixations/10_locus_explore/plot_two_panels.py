import sqlite3
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

matplotlib.rcParams.update({'font.size': 22})

TWON = 1e4
TIME_OFFSET = 5e4
ALPHA = 0.6
LINEWIDTH = 2


datasets = []
datasets_notfixed = []
for i in ['lowmu.22371142.sqlite3', 'midmu.1231116.sqlite3', 'himu.30884989.sqlite3']:
    with sqlite3.connect(i) as conn:
        df = pd.read_sql("select * from data", conn)
        datasets.append(df)
        df = pd.read_sql("select * from data_notfixed", conn)
        datasets_notfixed.append(df)

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(3, 3)

freqaxes = []
distanceaxes = []
saxes = []
saxes_inset = []

freqaxes.append(plt.subplot(gs[0, 0]))
freqaxes.append(plt.subplot(gs[0, 1], sharex=freqaxes[0], sharey=freqaxes[0]))
freqaxes.append(plt.subplot(gs[0, 2], sharex=freqaxes[0], sharey=freqaxes[0]))

distanceaxes.append(plt.subplot(gs[1, 0]))
distanceaxes.append(plt.subplot(
    gs[1, 1], sharex=distanceaxes[0], sharey=distanceaxes[0]))
distanceaxes.append(plt.subplot(
    gs[1, 2], sharex=distanceaxes[0], sharey=distanceaxes[0]))

saxes.append(plt.subplot(gs[2, 0]))
saxes.append(plt.subplot(gs[2, 1], sharex=saxes[0], sharey=saxes[0]))
saxes.append(plt.subplot(gs[2, 2], sharex=saxes[0], sharey=saxes[0]))

for i in range(len(saxes)):
    saxes_inset.append(inset_axes(saxes[i],
                                  width="50%",
                                  height=0.5,
                                  loc='upper right'))


tadapted_per_dataset = []

for i in range(len(datasets)):
    dset = datasets[i][5000*datasets[i].esize ** 2 >= 100]
    grps = dset.groupby('position')

    # Get the time when the replicate first reaches 0.9*z_o and 0.975*z_o
    first_group = list(grps)[0]
    idx = np.where(first_group[1].Gbar >= 0.9)[0]
    tadapted = first_group[1].generation.iloc[idx[0]] - TIME_OFFSET
    tadapted_per_dataset.append(tadapted)
    for g in grps:
        ftime = int(g[1].ftime.unique()[0] - TIME_OFFSET)
        origin = int(g[1].origin.unique()[0] - TIME_OFFSET)
        esize = g[1].esize.unique()[0]
        ls = '-'
        if esize < 0:
            ls = ':'
        if 5e3*esize*esize >= 100:
            xdata = g[1].generation-TIME_OFFSET
            Ng2 = 0.5*TWON*esize*esize
            label = r'$N\gamma^2$'+f"={Ng2:.3},\n{origin}, {ftime}"
            label = r'$N\gamma^2$'+f"={int(Ng2)},\n({origin}, {ftime})"
            if 5e3*esize*esize < 10:
                label = '_nolegend_'
            freqaxes[i].plot(xdata, g[1].nsamples/TWON, label=label,
                             linewidth=LINEWIDTH, ls=ls, alpha=ALPHA)

            distance = (g[1].G-1.0)  # /np.sqrt(g[1].VG)
            distance = (g[1].G-g[1].Gbar) / np.sqrt(g[1].VG)
            distanceaxes[i].plot(
                xdata, distance, linewidth=LINEWIDTH, ls=ls, alpha=ALPHA)
            s = (g[1].W - g[1].Wbar)/g[1].Wbar  # VS = 1
            saxes[i].plot(xdata, s, linewidth=LINEWIDTH, ls=ls, alpha=ALPHA)
            saxes_inset[i].plot(
                xdata, s, linewidth=LINEWIDTH, ls=ls, alpha=ALPHA)

    freqaxes[i].axvline(tadapted, color="grey",
                        ls="dashed", linewidth=LINEWIDTH, alpha=0.3)
    distanceaxes[i].axvline(tadapted, color="grey",
                            ls="dashed", linewidth=LINEWIDTH, alpha=0.3)
    saxes[i].axvline(tadapted, color="grey",
                     ls="dashed", linewidth=LINEWIDTH, alpha=0.3)
    saxes[i].axhline(2./TWON, color="grey",
                     ls="dotted", linewidth=LINEWIDTH, alpha=0.3)
    saxes_inset[i].axvline(tadapted, color="grey",
                           ls="dashed", linewidth=LINEWIDTH, alpha=0.3)
    saxes_inset[i].axhline(2./TWON, color="grey",
                           ls="dotted", linewidth=LINEWIDTH, alpha=0.3)

for axes in freqaxes[:1]:
    axes.legend(loc='best', prop={'size': 12},
                frameon=False)
for axes in freqaxes[1:]:
    axes.legend(loc='best', prop={'size': 12},
                frameon=False)

for axes in saxes_inset:
    # axes.set_xlim(saxes[0].get_xlim()[0],75)
    axes.set_xlim(20, 125)
    axes.set_ylim(0.5*saxes[0].get_ylim()[0], 500./TWON)
    for item in axes.get_xticklabels() + axes.get_yticklabels():
        item.set_fontsize(11)
    axes.xaxis.set_tick_params(width=0.25)
    axes.yaxis.set_tick_params(width=0.25)
    for loc in ['top', 'bottom', 'left', 'right']:
        axes.spines[loc].set_linewidth(0.25)

# Hide some axes!!!
for axes in freqaxes + distanceaxes:
    axes.get_xaxis().set_visible(False)

for axes in freqaxes[1:] + distanceaxes[1:] + saxes[1:]:
    axes.get_yaxis().set_visible(False)

saxes[1].set_xlabel("Generations since optimum shift")
saxes[0].set_ylabel(r'$\frac{\bar{w}_{a-} - \bar{w}}{\bar{w}}$')
# distanceaxes[0].set_ylabel(r'$\bar{z}_{a-} - z_o$')
distanceaxes[0].set_ylabel(r'$\frac{\bar{z}_{a-} - \bar{z}}{\sigma_z}$')
freqaxes[0].set_ylabel("Frequency")

distanceaxes[0].set_ylim(distanceaxes[0].get_ylim()[0], 12)
saxes[0].set_ylim(saxes[0].get_ylim()[0], 0.5)

for ax in saxes:
    ax.set_xticks([0, 100, 200])

mu = [r'$2.5\times 10^{-4}$', '0.001', '0.005']
for ax, m in zip(freqaxes, mu):
    ax.set_title(r'$\mu = $' + m, fontdict={'fontsize': 19})

plt.savefig("StrengthSelectionFixations.pdf")
plt.clf()

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(3, 3)


# Axes for non-fixations
freqaxes_notfixed = []
distanceaxes_notfixed = []
saxes_notfixed = []
freqaxes_notfixed.append(plt.subplot(gs[0, 0], sharey=freqaxes[0]))
freqaxes_notfixed.append(plt.subplot(
    gs[0, 1], sharex=freqaxes[0], sharey=freqaxes[0]))
freqaxes_notfixed.append(plt.subplot(
    gs[0, 2], sharex=freqaxes[0], sharey=freqaxes[0]))

distanceaxes_notfixed.append(plt.subplot(gs[1, 0], sharey=distanceaxes[0]))
distanceaxes_notfixed.append(plt.subplot(
    gs[1, 1], sharex=distanceaxes[0], sharey=distanceaxes[0]))
distanceaxes_notfixed.append(plt.subplot(
    gs[1, 2], sharex=distanceaxes[0], sharey=distanceaxes[0]))

saxes_notfixed.append(plt.subplot(gs[2, 0], sharey=saxes[0]))
saxes_notfixed.append(plt.subplot(gs[2, 1], sharex=saxes[0], sharey=saxes[0]))
saxes_notfixed.append(plt.subplot(gs[2, 2], sharex=saxes[0], sharey=saxes[0]))
# Plot notfixed
freq_thresh = 0.1*5000  # 5%
for i in range(len(datasets_notfixed)):
    tadapted = tadapted_per_dataset[i]
    print(tadapted)
    dset = datasets_notfixed[i][datasets_notfixed[i].origin <
                                tadapted + TIME_OFFSET]
    print(type(dset))
    print(len(dset.index), len(datasets_notfixed[i].index))
    grps = dset.groupby('position')
    # Get the time when the replicate first reaches 0.9*z_o and 0.975*z_o
    first_group = list(grps)[0]
    for g in grps:
        if g[1].nsamples.max() < freq_thresh:
            continue
        origin = int(g[1].origin.unique()[0] - TIME_OFFSET)
        esize = g[1].esize.unique()[0]
        ls = '-'
        if esize < 0:
            ls = ':'
        if 5e3*esize*esize >= 100:
            xdata = g[1].generation-TIME_OFFSET
            Ng2 = 0.5*TWON*esize*esize
            label = r'$N\gamma^2$'+f"={Ng2:.1}, {origin}, NA"
            label = '_nolegend_'
            freqaxes_notfixed[i].plot(xdata, g[1].nsamples/TWON, label=label,
                                      linewidth=LINEWIDTH, ls=ls, alpha=ALPHA)

            distance = (g[1].G-1.0)  # /np.sqrt(g[1].VG)
            distance = (g[1].G-g[1].Gbar) / np.sqrt(g[1].VG)
            distanceaxes_notfixed[i].plot(
                xdata, distance, linewidth=LINEWIDTH, ls=ls, alpha=ALPHA)
            s = (g[1].W - g[1].Wbar)/g[1].Wbar  # VS = 1
            saxes_notfixed[i].plot(
                xdata, s, linewidth=LINEWIDTH, ls=ls, alpha=ALPHA)
            # saxes_inset[i].plot(
            #     xdata, s, linewidth=LINEWIDTH, ls=ls, alpha=ALPHA)

    freqaxes_notfixed[i].axvline(tadapted, color="grey",
                                 ls="dashed", linewidth=LINEWIDTH, alpha=0.3)
    distanceaxes_notfixed[i].axvline(tadapted, color="grey",
                                     ls="dashed", linewidth=LINEWIDTH, alpha=0.3)
    saxes_notfixed[i].axvline(tadapted, color="grey", alpha=0.3,
                              ls="dashed", linewidth=LINEWIDTH)
    saxes_notfixed[i].axhline(2./TWON, color="grey", alpha=0.3,
                              ls="dotted", linewidth=LINEWIDTH)
    # saxes_inset[i].axvline(tadapted, color="grey",
    #                        ls="dashed", linewidth=LINEWIDTH)
    # saxes_inset[i].axhline(2./TWON, color="grey",
    #                        ls="dotted", linewidth=LINEWIDTH)


fig.align_ylabels([i[0] for i in [freqaxes, distanceaxes, saxes]])

for ax in freqaxes_notfixed[1:] + distanceaxes_notfixed[1:]:
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)

for ax in [freqaxes_notfixed[0], distanceaxes_notfixed[0]]:
    ax.get_xaxis().set_visible(False)

for ax in saxes_notfixed[1:]:
    ax.get_yaxis().set_visible(False)

for ax in saxes_notfixed:
    ax.set_xticks([0, 100, 200])

saxes_notfixed[1].set_xlabel("Generations since optimum shift")
saxes_notfixed[0].set_ylabel(r'$\frac{\bar{w}_{a-} - \bar{w}}{\bar{w}}$')
# distanceaxes[0].set_ylabel(r'$\bar{z}_{a-} - z_o$')
distanceaxes_notfixed[0].set_ylabel(
    r'$\frac{\bar{z}_{a-} - \bar{z}}{\sigma_z}$')
freqaxes_notfixed[0].set_ylabel("Frequency")

# Titles
for ax, m in zip(freqaxes_notfixed, mu):
    ax.set_title(r'$\mu = $' + m, fontdict={'fontsize': 19})

plt.savefig("StrengthSelectionLosses.pdf")
