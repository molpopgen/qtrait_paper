import sqlite3
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

TWON = 1e4
TIME_OFFSET = 5e4


datasets = []
for i in glob.glob('*.sqlite3'):
    with sqlite3.connect(i) as conn:
        df = pd.read_sql("select * from data", conn)
        datasets.append(df)

fig = plt.figure()
gs = gridspec.GridSpec(3, 3)

freqaxes = []
distanceaxes = []
saxes = []

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

for i in range(len(datasets)):
    grps = datasets[i].groupby('position')
    # Get the time when the replicate first reaches 0.9*z_o and 0.975*z_o
    first_group = list(grps)[0]
    idx = np.where(first_group[1].Gbar >= 0.9)[0]
    tadapted = idx[0]
    for g in grps:
        xdata = g[1].generation-TIME_OFFSET
        freqaxes[i].plot(xdata, g[1].nsamples/TWON)
        distance = (g[1].G-1.0)  # /np.sqrt(g[1].VG)
        distanceaxes[i].plot(xdata, distance)
        s = ((g[1].G - 1.0)**2)/2  # VS = 1
        esize = g[1].esize.unique()[0]
        ftime = int(g[1].ftime.unique()[0] - TIME_OFFSET)
        origin = int(g[1].origin.unique()[0] - TIME_OFFSET)
        saxes[i].plot(xdata, s, label=r'$\gamma$'+f"={esize:.3}, {origin}, {ftime}")

        freqaxes[i].axvline(tadapted, color="black", ls="dashed")
        distanceaxes[i].axvline(tadapted, color="black", ls="dashed")
        saxes[i].axvline(tadapted, color="black", ls="dashed")
        saxes[i].axhline(2./TWON, color="black", ls="dotted")

for axes in saxes:
    axes.legend(loc='upper right', prop={'size': 6})

# Hide some axes!!!
for axes in freqaxes + distanceaxes:
    axes.get_xaxis().set_visible(False)

for axes in freqaxes[1:] + distanceaxes[1:] + saxes[1:]:
    axes.get_yaxis().set_visible(False)

saxes[1].set_xlabel("Generations since optimum shift")
saxes[0].set_ylabel("Mean " + r'$s_{a-}$')
distanceaxes[0].set_ylabel("Mean distance")
freqaxes[0].set_ylabel("Frequency")

distanceaxes[0].set_ylim(-1, 0.3)
saxes[0].set_ylim(saxes[0].get_ylim()[0], 0.4)

fig.align_ylabels()
plt.savefig("StrengthOfSelection.pdf")
