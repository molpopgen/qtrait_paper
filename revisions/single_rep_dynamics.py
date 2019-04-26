import fwdpy11
import gzip
import sqlite3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np

NLOCI = 10
MAXGEN = 100


def get_gvalue_per_locus(pop):
    locus_boundaries = [(i, i + 11) for i in range(0, NLOCI * 11, 11)]
    amd = np.array(pop.ancient_sample_metadata, copy=False)
    nodes = np.array(pop.tables.nodes, copy=False)
    nt = nodes['time'][amd['nodes'][:, 0]]
    w = np.where((nt >= 5e4) & (nt <= 5e4 + MAXGEN))[0]
    unt = np.unique(nt[w])
    gv_per_locus = np.zeros(len(unt) * NLOCI).reshape(len(unt), NLOCI)
    dummy = 0
    for t in unt:
        wt = np.where(nt == t)[0]
        samples = amd['nodes'][wt].flatten()
        gv_per_locus_t = np.zeros(len(wt) * NLOCI).reshape(len(wt), NLOCI)
        vi = fwdpy11.ts.VariantIterator(pop.tables, pop.mutations, samples)
        for v in vi:
            r = v.record
            if pop.mutations[r.key].neutral is False:
                g = v.genotypes
                dg = g.reshape(len(wt), 2).sum(axis=1)
                pos = pop.mutations[r.key].pos
                for i, j in enumerate(locus_boundaries):
                    if pos < j[1]:
                        gv_per_locus_t[:, i] += (dg * pop.mutations[r.key].s)
                        break
        gv_per_locus[dummy, :] += gv_per_locus_t.mean(axis=0)
        dummy += 1
    return gv_per_locus


def get_data(case):
    with gzip.open("main_results/{}/replicate1.gz".format(case), "rb") as f:
        pop = fwdpy11.SlocusPop.load_from_pickle_file(f)

    gv_per_locus = get_gvalue_per_locus(pop)

    mg = 5e4 + MAXGEN
    with sqlite3.connect("main_results_{}_gscan.sqlite3".format(case)) as f:
        gscan_df = pd.read_sql(
            'select * from data where repid == 1 and window==5 and generation <= {}'.format(mg), f)

    with sqlite3.connect("main_results_{}_qtrait.sqlite3".format(case)) as f:
        qtrait_df = pd.read_sql(
            'select * from data where repid == 1 and generation <= {}'.format(mg), f)
    print(len(qtrait_df.index))

    return gscan_df, qtrait_df, gv_per_locus


def plot(traits, D, H, G, gscan_df, qtrait_df, gv_per_locus):
    traits.plot(qtrait_df.generation - 5e4, qtrait_df.zbar)
    traits.set_xlim(0, 100)

    gr = gscan_df.groupby(['locus'])
    for n, g in gr:
        D.plot(g.generation - 5e4, g.tajd, label=n, alpha=0.5)
        H.plot(g.generation - 5e4, g.hprime, label=n, alpha=0.5)
        G.plot(g.generation - 5e4, gv_per_locus[:, n], label=n, alpha=0.5)

    # D.legend(loc='best', fontsize='x-small')
    # H.legend(loc='best', fontsize='x-small')


fig = plt.figure()
gs = gridspec.GridSpec(4, 3)
lotraits = plt.subplot(gs[0, 0])
loG = plt.subplot(gs[1, 0], sharex=lotraits)
loD = plt.subplot(gs[2, 0], sharex=lotraits)
loH = plt.subplot(gs[3, 0], sharex=lotraits)
midtraits = plt.subplot(gs[0, 1])
midG = plt.subplot(gs[1, 1], sharex=midtraits, sharey=loG)
midD = plt.subplot(gs[2, 1], sharex=midtraits, sharey=loD)
midH = plt.subplot(gs[3, 1], sharex=midtraits, sharey=loH)
hitraits = plt.subplot(gs[0, 2])
hiG = plt.subplot(gs[1, 2], sharex=hitraits, sharey=loG)
hiD = plt.subplot(gs[2, 2], sharex=hitraits, sharey=loD)
hiH = plt.subplot(gs[3, 2], sharex=hitraits, sharey=loH)

gscan, traits, gv_per_locus = get_data('lomu')
plot(lotraits, loD, loH, loG, gscan, traits, gv_per_locus)

gscan, traits, gv_per_locus = get_data('midmu')
plot(midtraits, midD, midH, midG, gscan, traits, gv_per_locus)

gscan, traits, gv_per_locus = get_data('himu')
plot(hitraits, hiD, hiH, hiG, gscan, traits, gv_per_locus)

for axes in [lotraits, loD, loG, midtraits, midD, midG, hitraits, hiD, hiG]:
    axes.get_xaxis().set_visible(False)

for axes in [midtraits, midG, midD, midH, hitraits, hiG, hiD, hiH]:
    axes.get_yaxis().set_visible(False)

lotraits.set_ylabel(r'$\bar{z}$')
lotraits.set_ylim(0, 1)
loD.set_ylabel("Tajima's D")
loH.set_ylabel("H'")
loG.set_ylabel(r'$\bar{z}$ per locus')
loG.set_ylim(0,1)
midH.set_xlabel("Generations since optimum shift")

lotraits.set_title(r'$\mu = 2.5\times 10^{-4}$')
midtraits.set_title(r'$\mu = 1\times 10^{-3}$')
hitraits.set_title(r'$\mu = 5\times 10^{-3}$')

fig.align_ylabels()
plt.savefig('SingleRepGenomeScan.pdf')
