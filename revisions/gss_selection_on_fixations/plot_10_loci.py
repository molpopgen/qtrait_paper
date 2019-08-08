import sqlite3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.gridspec as gridspec

matplotlib.rcParams.update({'font.size': 22})

NLOCI = 10
LOCUS_LENGTH = 11

LOCUS_BOUNDARIES = [(i, i + 11) for i in range(0, 10 * 11, 11)]
SREGIONS = [(i[0]+5, i[0]+6) for i in LOCUS_BOUNDARIES]

POPSIZE = 1000

fig = plt.figure(figsize=(15, 15))
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

files = ['pop10loci_100.sqlite3',
         'pop10loci_1000.sqlite3', 'pop10loci_10000.sqlite3']


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
    colors = ['w']*len(Ns)
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
        fixations = pd.read_sql(
            'select * from fixations where origin <= 10100 and fixation > 10000', conn)

    fixations['alpha'] = esize_to_alpha(fixations)
    fixations['color'] = esize_to_color(fixations)
    grps = data.groupby(['threshold'])
    for n, g in grps:
        Daxes[i].plot(g.left_edge, g.D, alpha=0.5)
        Haxes[i].plot(g.left_edge, g.Hp, alpha=0.5)
        Hdivaxes[i].plot(g.left_edge, g.hdiv, alpha=0.5, label=n)

        for s in SREGIONS:
            Daxes[i].axvspan(*s, facecolor='lightgray', alpha=0.1)
            Haxes[i].axvspan(*s, facecolor='lightgray', alpha=0.1)
            Hdivaxes[i].axvspan(*s, facecolor='lightgray', alpha=0.1)

for i, f in enumerate(files):
    with sqlite3.connect(f) as conn:
        fixations = pd.read_sql(
            'select * from fixations where origin <= 10100 and fixation > 10000', conn)

    fixations['alpha'] = esize_to_alpha(fixations)
    fixations['color'] = esize_to_color(fixations)

    for p, a, c in zip(fixations.position, fixations.alpha, fixations.color):
        lim = Daxes[0].get_ylim()
        Daxes[i].plot(p, Daxes[0].dataLim.y1, color=c, marker='v')
        lim = Haxes[0].get_ylim()
        Haxes[i].plot(p, Haxes[0].dataLim.y1, color=c, marker='v')
        lim = Hdivaxes[0].get_ylim()
        Hdivaxes[i].plot(p, 1.0, color=c, marker='v')

for ax in Daxes[1:] + Haxes[1:]:
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)

for ax in Hdivaxes[1:]:
    ax.get_yaxis().set_visible(False)

for ax in Hdivaxes:
    ax.legend(loc='best', frameon=False)

Daxes[0].set_ylabel("Tajima's D")
Haxes[0].set_ylabel("H'")
Hdivaxes[0].set_ylabel("Haplotype diversity")

plt.tight_layout()

plt.savefig('test.pdf')
