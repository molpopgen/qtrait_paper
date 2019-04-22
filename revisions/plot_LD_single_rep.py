import sqlite3
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np

THRESHOLD = np.sqrt(50. / 5e3)

fig = plt.figure()
gs = gridspec.GridSpec(3, 3)
freqaxes = []
ldaxes = []

freqaxes.append(fig.add_subplot(gs[0, 0]))
freqaxes.append(fig.add_subplot(gs[0, 1],sharex=freqaxes[0]))
freqaxes.append(fig.add_subplot(gs[0, 2],sharex=freqaxes[0]))

ldaxes.append(fig.add_subplot(gs[1, 0], sharex=freqaxes[0]))
ldaxes.append(fig.add_subplot(gs[2, 0], sharex=freqaxes[0], sharey=ldaxes[0]))

ldaxes.append(fig.add_subplot(gs[1, 1], sharex=freqaxes[1], sharey=ldaxes[0]))
ldaxes.append(fig.add_subplot(gs[2, 1], sharex=freqaxes[1], sharey=ldaxes[0]))

ldaxes.append(fig.add_subplot(gs[1, 2], sharex=freqaxes[2], sharey=ldaxes[0]))
ldaxes.append(fig.add_subplot(gs[2, 2], sharex=freqaxes[2], sharey=ldaxes[0]))

directories = ['lowrec', 'main_results', 'hirec']
replicates=[199]*3
replicates=[189]*3
replicates=[169]*3
replicates=[1]*3
for i, d in enumerate(directories):
    conn = sqlite3.connect('{}/himu/replicate{}_ld.sqlite3'.format(d,replicates[i]))
    df = pd.read_sql(
        'select * from data where abs(pos2-pos1) < 1.0 and generation < 60000 and abs(e1) >= {}'.format(THRESHOLD), conn)

    # plot freqs in first plot
    groups = df.groupby(['pos1'])
    for n, g in groups:
        print(g.e1.unique())
        freqaxes[i].plot(g.generation - 5e4, g.c1 / 1e4,
                         label='{:0.3f}'.format(g.e1.unique()[0]))
        freqaxes[i].legend(loc='lower right', fontsize='small', edgecolor=None,
                framealpha=0.0)

    dplotneg = ldaxes[2 * i]
    dplotpos = ldaxes[2 * i + 1]
    minc2 = 0.10 * 1e4
    alpha = 1
    for n, g in groups:
        g2 = g[(g.c2 >= minc2) & (g.c2 < 1e4 - minc2) & (g.e2 < 0)]
        g2 = g2.groupby(['generation']).mean().reset_index()
        dplotneg.scatter(g2.generation - 5e4, g2.D, s=0.05, alpha=alpha)
        g2 = g[(g.c2 >= minc2) & (g.c2 < 1e4 - minc2) & (g.e2 > 0)]
        g2 = g2.groupby(['generation']).mean().reset_index()
        dplotpos.scatter(g2.generation - 5e4, g2.D, s=0.05, alpha=alpha)

# hide axes
for i in freqaxes:
    i.axes.get_xaxis().set_visible(False)
for i in freqaxes[1:]:
    i.axes.get_yaxis().set_visible(False)
for i in ldaxes[0:5:2]:
    i.axes.get_xaxis().set_visible(False)
for i in ldaxes[2:]:
    i.axes.get_yaxis().set_visible(False)

# axis labels
ldaxes[3].set_xlabel("Generations since optimum shift")
freqaxes[0].set_ylabel("Frequency")
ldaxes[0].set_ylabel(r'$\bar{D}, \gamma < 0$')
ldaxes[1].set_ylabel(r'$\bar{D}, \gamma > 0$')

# column titles
for i,j in zip(freqaxes, [1e2,1e3,1e4]):
    i.set_title(r'$\rho = {}$'.format(j))

fig.align_ylabels()
plt.savefig('LDsinglerep.pdf')
