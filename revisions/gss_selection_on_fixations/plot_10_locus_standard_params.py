import sqlite3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.gridspec as gridspec

matplotlib.rcParams.update({'font.size': 22})

CMAP = matplotlib.cm.get_cmap('viridis')
LINECOLORS = [j for j in reversed([CMAP(1. - i / 3) for i in range(1, 4)])]
LINECOLORS = [CMAP(1. - i / 3) for i in range(1, 4)]

NLOCI = 10
LOCUS_LENGTH = 11

LOCUS_BOUNDARIES = [(i, i + 11) for i in range(0, 10 * 11, 11)]
SREGIONS = [(i[0] + 5, i[0] + 6) for i in LOCUS_BOUNDARIES]

POPSIZE = 5000

fig = plt.figure(figsize=(25, 15))
gs = gridspec.GridSpec(3, 3)
D = plt.subplot(gs[0, 0])
D1 = plt.subplot(gs[0, 1], sharey=D)
D2 = plt.subplot(gs[0, 2], sharey=D)

H = plt.subplot(gs[1, 0])
H1 = plt.subplot(gs[1, 1], sharey=H)
H2 = plt.subplot(gs[1, 2], sharey=H)

Hdiv = plt.subplot(gs[2, 0])
Hdiv1 = plt.subplot(gs[2, 1], sharey=Hdiv)
Hdiv2 = plt.subplot(gs[2, 2], sharey=Hdiv)

Daxes = [D, D1, D2]
Haxes = [H, H1, H2]
Hdivaxes = [Hdiv, Hdiv1, Hdiv2]

files = ['pop_lomu_1000.sqlite3', 'pop_midmu_1000.sqlite3',
         'pop_himu_1000.sqlite3']


def esize_to_alpha(fixations):
    Ns = POPSIZE * fixations.esize**2
    logrange = np.logspace(0, 1, num=4)
    logrange /= logrange[-1]
    alphas = np.zeros(len(Ns))
    for i, j in enumerate(Ns):
        if j <= 1:
            alphas[i] = logrange[0]
        elif j <= 10:
            alphas[i] = logrange[1]
        elif j <= 100:
            alphas[i] = logrange[2]
        else:
            alphas[i] = logrange[3]
    return alphas


def esize_to_color(fixations):
    Ns = POPSIZE * fixations.esize**2
    colors = ['w'] * len(Ns)
    for i, j in enumerate(Ns):
        if j <= 1:
            colors[i] = 'w'
        elif j <= 10:
            colors[i] = 'y'
        elif j <= 100:
            colors[i] = 'c'
        else:
            colors[i] = 'm'
    return colors


for i, f in enumerate(files):
    with sqlite3.connect(f) as conn:
        data = pd.read_sql('select * from data', conn)

    grps = data.groupby(['threshold'])
    c = 1
    for n, g in grps:
        lcolor = CMAP(1. - c / 3.)
        lcolor = CMAP(c / 3.)
        lcolor = LINECOLORS[c - 1]
        lalpha = 1.2 - c / 5.
        loc = np.unique(g.locus)
        for lli, ll in enumerate(loc):
            idx = np.where(g.locus == ll)[0]
            le = np.array(g.left_edge)
            pi = np.array(g.pi)
            Hp = np.array(g.Hp)
            hd = np.array(g.hdiv)
            Daxes[i].plot(le[idx], pi[idx],
                          alpha=lalpha, color=lcolor, linewidth=1)
            Haxes[i].plot(le[idx], Hp[idx],
                          alpha=lalpha, color=lcolor, linewidth=1)
            if lli == 0:
                Hdivaxes[i].plot(le[idx], hd[idx], alpha=lalpha,
                                 color=lcolor, label=n, linewidth=1)
            else:
                Hdivaxes[i].plot(le[idx], hd[idx], alpha=lalpha,
                                 color=lcolor, label='_nolegend_', linewidth=1)
        c += 1

        for s in SREGIONS:
            Daxes[i].axvspan(*s, facecolor='lightgray', alpha=0.1)
            Haxes[i].axvspan(*s, facecolor='lightgray', alpha=0.1)
            Hdivaxes[i].axvspan(*s, facecolor='lightgray', alpha=0.1)

for i, f in enumerate(files):
    with sqlite3.connect(f) as conn:
        fixations = pd.read_sql(
            'select * from fixations where origin <= 50100 and fixation > 50000', conn)

    fixations['alpha'] = esize_to_alpha(fixations)
    fixations['color'] = esize_to_color(fixations)

    for p, a, c, o in zip(fixations.position, fixations.alpha, fixations.color, fixations.origin):
        marker = 'v'
        if o < 10 * POPSIZE:
            marker = '^'
        lim = Daxes[0].get_ylim()
        Daxes[i].plot(p, Daxes[0].dataLim.y1, color=c, marker=marker, markersize=12)
        lim = Haxes[0].get_ylim()
        Haxes[i].plot(p, Haxes[0].dataLim.y1, color=c, marker=marker, markersize=12)
        lim = Hdivaxes[0].get_ylim()
        Hdivaxes[i].plot(p, 1.0, color=c, marker=marker, markersize=12)

for i, j in zip(Daxes, [2.5e-4, 1e-3, 5e-3]):
    i.set_title(r'$\mu =$' + '{:.2}'.format(j))

for ax in Daxes + Haxes:
    ax.get_xaxis().set_visible(False)

for ax in Daxes[1:] + Haxes[1:]:
    ax.get_yaxis().set_visible(False)

for ax in Hdivaxes[1:]:
    ax.get_yaxis().set_visible(False)

for ax in Hdivaxes:
    ax.legend(loc='best', frameon=False)

Daxes[0].set_ylabel(r'$\pi$')
Haxes[0].set_ylabel("H'")
Hdivaxes[0].set_ylabel("Haplotype diversity")

for ax in Hdivaxes:
    ax.set_xticks([i[0] + 5.5 for i in LOCUS_BOUNDARIES])
    ax.set_xticklabels([i + 1 for i in range(len(LOCUS_BOUNDARIES))])

Hdivaxes[1].set_xlabel('Locus mid-point')

plt.tight_layout()

plt.savefig('SpatialTemporalVariation.pdf')
