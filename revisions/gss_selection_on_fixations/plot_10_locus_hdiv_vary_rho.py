import sqlite3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.lines

matplotlib.rcParams.update({'font.size': 22})

CMAP = matplotlib.cm.get_cmap('viridis')
LINECOLORS = [j for j in reversed([CMAP(1. - i / 3) for i in range(1, 4)])]
LINECOLORS = [CMAP(1. - i / 3) for i in range(1, 4)]

NLOCI = 10
LOCUS_LENGTH = 11

LOCUS_BOUNDARIES = [(i, i + 11) for i in range(0, 10 * 11, 11)]
SREGIONS = [(i[0] + 5, i[0] + 6) for i in LOCUS_BOUNDARIES]

POPSIZE = 5000


filessets = [
    ['pop_lomu_100.sqlite3', 'pop_midmu_100.sqlite3',
        'pop_himu_100.sqlite3'],
    ['pop_lomu_1000.sqlite3', 'pop_midmu_1000.sqlite3',
        'pop_himu_1000.sqlite3'],
    ['pop_lomu_10000.sqlite3', 'pop_midmu_10000.sqlite3',
     'pop_himu_10000.sqlite3']
]
filessets = [
    ['pop_lomu_100.sqlite3', 'pop_lomu_1000.sqlite3',
     'pop_lomu_10000.sqlite3'],
    ['pop_midmu_100.sqlite3', 'pop_midmu_1000.sqlite3',
     'pop_midmu_10000.sqlite3'],
    ['pop_himu_100.sqlite3', 'pop_himu_1000.sqlite3',
     'pop_himu_10000.sqlite3']
]

mutrates = [2.5e-4, 1e-3, 5e-3]


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


fig = plt.figure(figsize=(25, 15))
gs = gridspec.GridSpec(3, 3)
D = plt.subplot(gs[0, 0])
D1 = plt.subplot(gs[0, 1], sharey=D)
D2 = plt.subplot(gs[0, 2], sharey=D)

H = plt.subplot(gs[1, 0], sharey=D)
H1 = plt.subplot(gs[1, 1], sharey=D)
H2 = plt.subplot(gs[1, 2], sharey=D)

Hdiv = plt.subplot(gs[2, 0], sharey=D)
Hdiv1 = plt.subplot(gs[2, 1], sharey=D)
Hdiv2 = plt.subplot(gs[2, 2], sharey=D)

Daxes = [D, D1, D2]
Haxes = [H, H1, H2]
Hdivaxes = [Hdiv, Hdiv1, Hdiv2]
axes = [Daxes, Haxes, Hdivaxes]


for fset, mu in zip(enumerate(filessets), mutrates):
    fsetindex, files = fset
    for i, f in enumerate(files):
        print(f)
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
                hd = np.array(g.hdiv)
                label = n
                if lli > 0:
                    label = '_nolegend_'
                axes[fsetindex][i].scatter(le[idx], hd[idx], label=label,
                                           color=lcolor, marker='.')
                # Haxes[i].plot(g.left_edge, g.Hp, alpha=lalpha,
                #               color=lcolor, linewidth=1)
                # Hdivaxes[i].plot(g.left_edge, g.hdiv, alpha=lalpha,
                #                  color=lcolor, label=n, linewidth=1)
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

        for p, a, c, o in zip(fixations.position, fixations.alpha,
                              fixations.color, fixations.origin):
            marker = 'v'
            if o <= 10 * POPSIZE:
                marker = '^'
            axes[fsetindex][i].plot(
                p, 1.05, color=c, marker=marker, markersize=12)

for i, j in zip(Daxes,  [1e2, 1e3, 1e4]):
    i.set_title(r'$\rho =$' + '{:.2}'.format(j))

for ax in Daxes + Haxes:
    ax.get_xaxis().set_visible(False)

for ax in Daxes[1:] + Haxes[1:]:
    ax.get_yaxis().set_visible(False)

for ax in Hdivaxes[1:]:
    ax.get_yaxis().set_visible(False)

thresholds = [0.1, 0.5, 0.9]
handles = [matplotlib.lines.Line2D([0], [0], marker='.', markersize=18, label=thresholds[i - 1],
                                   linestyle='',
                                   color=CMAP(1 - i / 3)) for i in range(1, 4)]
for a, m in zip(axes, mutrates):
    for ax in a:
        ax.text(-3, 0.025, r'$\mu = $' + '{:.2}'.format(m))
        ax.legend(loc='best', bbox_to_anchor=(
            0.0, 0.0, 0.25, 0.5), frameon=False,
            handles=handles)


# for ax in Hdivaxes:
#     ax.legend(loc='best', frameon=False)

Daxes[0].set_ylabel("Haplotype diversity")
Haxes[0].set_ylabel("Haplotype diversity")
Hdivaxes[0].set_ylabel("Haplotype diversity")

for ax in Hdivaxes:
    ax.set_xticks([i[0] + 5.5 for i in LOCUS_BOUNDARIES])
    ax.set_xticklabels([i + 1 for i in range(len(LOCUS_BOUNDARIES))])

Hdivaxes[1].set_xlabel('Locus mid-point')

plt.tight_layout()

plt.savefig('HapDivVaryRho.pdf')
