import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
from cycler import cycler
import database
from sys import exit
import numpy as np

dbname = "data/spectra3.db"

ETminRaw = 10
ETmaxRaw = 32.5
ETminRen = 10
ETmaxRen = 18

output = "pdf"
renlist = ("raw", "rentails", "renloc")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

params = {'legend.fontsize': 8, 'lines.markersize':2.5, 'lines.marker':"o"}
plt.rcParams.update(params)
plt.rc('axes', prop_cycle=(cycler('linestyle', ['-', '--', ':', '-.'])))

color = {1:"b", -1:"r"}

neigs = 3

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 1.5


def plotvsET(ETlistRaw, ETlistRen):

    ETlist = {"raw":ETlistRaw, "renloc":ETlistRaw, "rentails":ETlistRen}

    db = database.Database(dbname)

    spectrum = {k: {ren:
        np.array([db.getEigs(k,ren,g,L,ET) for ET in ETlist[ren]]).transpose()
        for ren in renlist} for k in (-1,1)}


    # VACUUM
    plt.figure(1)
    for ren in renlist:
        data = spectrum[1][ren][0]
        label = "ren="+ren
        plt.plot(ETlist[ren], data, label=label, color='k')

    # SPECTRUM
    plt.figure(2)
    for k, n0 in ((1,1),(-1,0)):

        for i in range(n0,neigs):

            for ren in renlist:
                data = spectrum[k][ren][i]-spectrum[1][ren][0]
                if i==n0:
                    label="k={}, ren={}".format(k,ren)
                else:
                    label = None
                plt.plot(ETlist[ren], data, label=label, color=color[k])

            plt.gca().set_prop_cycle(None)

    # MASS
    plt.figure(3)
    for ren in renlist:
        data = spectrum[-1][ren][0]-spectrum[1][ren][0]
        label="ren={}".format(ren)
        plt.plot(ETlist[ren], data, label=label, color="k")

    plt.gca().set_prop_cycle(None)


argv = sys.argv


# if len(argv) < 5:
    # print(argv[0], "<L> <g> <ETmin> <ETmax>")
    # sys.exit(-1)

# L = float(argv[1])
# g = float(argv[2])
# ETmin = float(argv[3])
# ETmax = float(argv[4])

if len(argv) < 3:
    print(argv[0], "<L> <g>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])

ETlistRaw = scipy.linspace(ETminRaw, ETmaxRaw, (ETmaxRaw-ETminRaw)*2+1)
ETlistRen = scipy.linspace(ETminRen, ETmaxRen, (ETmaxRen-ETminRen)*2+1)

print("g=", g)


plotvsET(ETlistRaw, ETlistRen)


title = r"$g$={0:.1f}, $L$={1:.1f}".format(g,L)
fname = "g={0:.1f}_L={1:.1f}.{2}".format(g,L,output)


# VACUUM
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_{T}$")
plt.ylabel(r"$\mathcal{E}_0$")
plt.xlim(min(min(ETlistRaw),min(ETlistRen)), max(max(ETlistRaw),max(ETlistRen)))
plt.legend(loc=1)
plt.savefig("E0vsET_"+fname)


# MASSES
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_{T}$")
plt.ylabel(r"$\mathcal{E}_I - \mathcal{E}_0$")
plt.xlim(min(min(ETlistRaw),min(ETlistRen)), max(max(ETlistRaw),max(ETlistRen)))
plt.legend(loc=2)
plt.savefig("specvsET_"+fname)

# MASS
plt.figure(3, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_{T}$")
plt.ylabel(r"$\mathcal{E}_1 - \mathcal{E}_0$")
plt.xlim(min(min(ETlistRaw),min(ETlistRen)), max(max(ETlistRaw),max(ETlistRen)))
plt.legend(loc=2)
plt.savefig("massvsET_"+fname)
