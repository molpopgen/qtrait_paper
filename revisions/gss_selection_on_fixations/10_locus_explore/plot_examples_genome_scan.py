import sys
import sqlite3
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

matplotlib.rcParams.update({'font.size': 22})

TWON = 1e4
TIME_OFFSET = 5e4
ALPHA = 0.8
LINEWIDTH = 2
FREQ_THRESH = 0.05 * TWON  # 5%
COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]
LOCUS_BOUNDARIES = [(i, i + 11) for i in range(0, 10 * 11, 11)]
SREGIONS = [(i[0] + 5, i[0] + 6) for i in LOCUS_BOUNDARIES]


def get_locus(df):
    pos = df.position.iloc[0]
    for i, j in enumerate(SREGIONS):
        if pos >= j[0] and pos < j[1]:
            return i
    raise RuntimeError("could not get locus ID")


datasets = []
datasets_notfixed = []
gscan = []
for i in ['lowmu.22371142.sqlite3', 'midmu.1231116.sqlite3', 'himu.30884989.sqlite3']:
    with sqlite3.connect(i) as conn:
        df = pd.read_sql("select * from data", conn)
        datasets.append(df)
        df = pd.read_sql("select * from data_notfixed", conn)
        datasets_notfixed.append(df)
    ii = i.replace("sqlite3", "genome_scan.sqlite3")
    with sqlite3.connect(ii) as conn:
        df = pd.read_sql("select * from data_sweep", conn)
        gscan.append(df)

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(5, 3)

fixaxes = []
lossaxes = []
Daxes = []
Hprimeaxes = []
Hdivaxes = []

fixaxes.append(plt.subplot(gs[0, 0]))
fixaxes.append(plt.subplot(gs[0, 1], sharex=fixaxes[0], sharey=fixaxes[0]))
fixaxes.append(plt.subplot(gs[0, 2], sharex=fixaxes[0], sharey=fixaxes[0]))

lossaxes.append(plt.subplot(gs[1, 0], sharey=fixaxes[0]))
lossaxes.append(plt.subplot(gs[1, 1], sharex=fixaxes[0], sharey=fixaxes[0]))
lossaxes.append(plt.subplot(gs[1, 2], sharex=fixaxes[0], sharey=fixaxes[0]))

Daxes.append(plt.subplot(gs[2, 0]))
Daxes.append(plt.subplot(gs[2, 1], sharex=Daxes[0], sharey=Daxes[0]))
Daxes.append(plt.subplot(gs[2, 2], sharex=Daxes[0], sharey=Daxes[0]))

Hprimeaxes.append(plt.subplot(gs[3, 0]))
Hprimeaxes.append(plt.subplot(
    gs[3, 1], sharex=Hprimeaxes[0], sharey=Hprimeaxes[0]))
Hprimeaxes.append(plt.subplot(
    gs[3, 2], sharex=Hprimeaxes[0], sharey=Hprimeaxes[0]))

Hdivaxes.append(plt.subplot(gs[4, 0]))
Hdivaxes.append(plt.subplot(gs[4, 1], sharex=Hdivaxes[0], sharey=Hdivaxes[0]))
Hdivaxes.append(plt.subplot(gs[4, 2], sharex=Hdivaxes[0], sharey=Hdivaxes[0]))

tadapted_per_dataset = []

for i in range(len(datasets)):
    grps = datasets[i].groupby('position')

    # Get the time when the replicate first reaches 0.9*z_o
    first_group = list(grps)[0]
    idx = np.where(first_group[1].Gbar >= 0.9)[0]
    tadapted = first_group[1].generation.iloc[idx[0]] - TIME_OFFSET
    tadapted_per_dataset.append(tadapted)
    loci_represented = []
    focal_locus = []  # Loci w/fixations with Ngamma^2 >= 100
    for g in grps:
        ftime = int(g[1].ftime.unique()[0] - TIME_OFFSET)
        origin = int(g[1].origin.unique()[0] - TIME_OFFSET)
        esize = g[1].esize.unique()[0]
        Ng2 = 0.5*TWON*esize*esize
        if Ng2 >= 1:
            locus = get_locus(g[1])
            loci_represented.append(locus)
            xdata = g[1].generation-TIME_OFFSET
            label = r'$N\gamma^2$'+f"={Ng2:.3},\n{origin}, {ftime}"
            label = r'$N\gamma^2$'+f"={int(Ng2)},\n({origin}, {ftime})"
            if 5e3*esize*esize < 10:
                label = '_nolegend_'
            alpha = ALPHA
            ls = '-'
            if Ng2 >= 100:
                focal_locus.append(locus)
            if Ng2 < 100 and Ng2 >= 10:
                ls = '--'
                alpha = ALPHA/2.
            elif Ng2 < 10:
                ls = ":"
                alpha = ALPHA/2.
            fixaxes[i].plot(xdata, g[1].nsamples/TWON, label=label,
                            color=COLORS[locus],
                            linewidth=LINEWIDTH, ls=ls, alpha=alpha)

    dset = datasets_notfixed[i][datasets_notfixed[i].origin <
                                tadapted + TIME_OFFSET]
    grps = dset.groupby('position')
    for g in grps:
        if g[1].nsamples.max() < FREQ_THRESH:
            continue
        origin = int(g[1].origin.unique()[0] - TIME_OFFSET)
        esize = g[1].esize.unique()[0]
        ls = '-'
        if esize < 0:
            ls = ':'
        if 5e3*esize*esize >= 1:
            locus = get_locus(g[1])
            loci_represented.append(locus)
            xdata = g[1].generation-TIME_OFFSET
            Ng2 = 0.5*TWON*esize*esize
            label = r'$N\gamma^2$'+f"={Ng2:.1}, {origin}, NA"
            label = '_nolegend_'
            alpha = ALPHA
            ls = '-'
            if Ng2 < 100 and Ng2 >= 10:
                ls = '--'
                alpha = ALPHA/2.
            elif Ng2 < 10:
                ls = ":"
                alpha = ALPHA/2.

            lossaxes[i].plot(xdata, g[1].nsamples/TWON, label=label, color=COLORS[locus],
                             linewidth=LINEWIDTH, ls=ls, alpha=alpha)

    dset = gscan[i]
    grps = dset.groupby('locus')
    for n, g in grps:
        if n in loci_represented:
            xdata = g.generation-TIME_OFFSET
            if n in focal_locus:
                alpha = ALPHA
                zorder = 10
            else:
                alpha = 0.25
                zorder = 1
            Daxes[i].plot(xdata, g.D, color=COLORS[n], linewidth=LINEWIDTH, alpha=alpha,
                          zorder=zorder)
            Hprimeaxes[i].plot(xdata, g.Hp, color=COLORS[n],
                               linewidth=LINEWIDTH, alpha=alpha,
                               zorder=zorder)
            Hdivaxes[i].plot(xdata, g.hdiv, color=COLORS[n],
                             linewidth=LINEWIDTH, alpha=alpha,
                             zorder=zorder)

    for ax in [fixaxes[i], lossaxes[i], Daxes[i], Hprimeaxes[i], Hdivaxes[i]]:
        ax.axvline(tadapted, color="grey",
                   ls="dashed", linewidth=LINEWIDTH, alpha=0.3)

for ax in fixaxes + lossaxes + Daxes + Hprimeaxes:
    ax.get_xaxis().set_visible(False)

for ax in fixaxes[1:] + lossaxes[1:] + Daxes[1:] + Hprimeaxes[1:] + Hdivaxes[1:]:
    ax.get_yaxis().set_visible(False)

for ax in Hdivaxes:
    ax.set_xticks([0, 100, 200])

fixaxes[0].set_ylabel("Frequency", fontdict={'fontsize': 16})
lossaxes[0].set_ylabel("Frequency", fontdict={'fontsize': 16})
Daxes[0].set_ylabel("Tajima's D", fontdict={'fontsize': 16})
Hprimeaxes[0].set_ylabel("H'", fontdict={'fontsize': 16})
Hdivaxes[0].set_ylabel("Haplotype\ndiversity", fontdict={'fontsize': 16})
Hdivaxes[1].set_xlabel("Generations since optimum shift",
                       fontdict={'fontsize': 16})

mu = [r'$2.5\times 10^{-4}$', '0.001', '0.005']
for ax, m in zip(fixaxes, mu):
    ax.set_title(r'$\mu = $' + m, fontdict={'fontsize': 19})

fig.align_ylabels([i[0]
                   for i in [fixaxes, lossaxes, Daxes, Hprimeaxes, Hdivaxes]])

plt.savefig("StrengthSelectionExampleGenomeScan.pdf")
