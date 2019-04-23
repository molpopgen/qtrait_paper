import sqlite3
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd


def get_data(case):
    with sqlite3.connect("main_results_{}_gscan.sqlite3".format(case)) as f:
        gscan_df = pd.read_sql(
            'select * from data where repid == 1 and window==5', f)

    with sqlite3.connect("main_results_{}_qtrait.sqlite3".format(case)) as f:
        qtrait_df = pd.read_sql('select * from data where repid == 1', f)

    return gscan_df, qtrait_df


def plot(traits, D, H, gscan_df, qtrait_df):
    traits.plot(qtrait_df.generation - 5e4, qtrait_df.zbar)
    traits.set_xlim(0, 100)

    gr = gscan_df.groupby(['locus'])
    for n, g in gr:
        D.plot(g.generation - 5e4, g.tajd, label=n, alpha=0.5)
        H.plot(g.generation - 5e4, g.hprime, label=n, alpha=0.5)

    # D.legend(loc='best', fontsize='x-small')
    # H.legend(loc='best', fontsize='x-small')


fig = plt.figure()
gs = gridspec.GridSpec(3, 3)
lotraits = plt.subplot(gs[0, 0])
loD = plt.subplot(gs[1, 0], sharex=lotraits)
loH = plt.subplot(gs[2, 0], sharex=lotraits)
midtraits = plt.subplot(gs[0, 1])
midD = plt.subplot(gs[1, 1], sharex=midtraits, sharey=loD)
midH = plt.subplot(gs[2, 1], sharex=midtraits, sharey=loH)
hitraits = plt.subplot(gs[0, 2])
hiD = plt.subplot(gs[1, 2], sharex=hitraits, sharey=loD)
hiH = plt.subplot(gs[2, 2], sharex=hitraits, sharey=loH)

gscan, traits = get_data('lomu')
plot(lotraits, loD, loH, gscan, traits)

gscan, traits = get_data('midmu')
plot(midtraits, midD, midH, gscan, traits)

gscan, traits = get_data('himu')
plot(hitraits, hiD, hiH, gscan, traits)

for axes in [lotraits, loD, midtraits, midD, hitraits, hiD]:
    axes.get_xaxis().set_visible(False)

for axes in [midtraits, midD, midH, hitraits, hiD, hiH]:
    axes.get_yaxis().set_visible(False)

lotraits.set_ylabel(r'$\bar{z}$')
loD.set_ylabel("Tajima's D")
loH.set_ylabel("H'")

lotraits.set_title(r'$\mu = 2.5\times 10^{-4}$')
midtraits.set_title(r'$\mu = 1\times 10^{-3}$')
hitraits.set_title(r'$\mu = 5\times 10^{-3}$')

fig.align_ylabels()
plt.savefig('SingleRepGenomeScan.pdf')
