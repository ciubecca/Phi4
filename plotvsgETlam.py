import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
from cycler import cycler
import database
from sys import exit, argv

output = "png"

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

g2 = -0.75

lammin = 4
ETmin = 10

klist = (1,-1)
neigs = 4

if len(argv) < 4:
    print("{} <L> <Emax> <Lambda>".format(argv[0]))
    exit(1)

L = float(argv[1])
Emax = float(argv[2])
Lambda = float(argv[3])


g4list = np.linspace(0,30,15)
lamlist = np.linspace(lammin, Lambda, 10)
ETlist = np.linspace(ETmin, Emax, 10)

marker = 'o'
markersize = 2.5

db = database.Database()


def plotvsET(L, lam, g2, g4, ETlist):

    spectrum = {k:[] for k in klist}
    masses = {}

    for k in klist:
        for ET in ETlist:

            approxQuery = {"g4":g4, g2:"g2", "L":L, "ET":ET, "Lambda":lam}
            exactQuery = {"k": k, "neigs":neigs}
            boundQuery = {}

            try:
                spectrum[k].append(db.getObjList('spec', exactQuery,
                    approxQuery, orderBy="date")[0])

            except IndexError:
                print("Not found:", exactQuery, approxQuery)
                exit(-1)

    for k in klist:
        spectrum[k] = array(spectrum[k])

    # Mass
    for k in (-1,1):
        if k==1:
            imin = 1
        else:
            imin = 0
        masses[k] = spectrum[k][:,imin]-spectrum[1][:,0]

    # SPECTRUM
    for k in (1):
        plt.figure(fignum(k))
        # for i in range(neigs):
        for i in range(1):
            data = spectrum[k][:,i]/L
            label = "lam={}, g={}".format(lam,g)
            plt.plot(ETlist, data, label=label, marker=marker, markersize=markersize)
        plt.gca().set_prop_cycle(None)

    # MASS
    plt.figure(2)
    # for k in (-1,1):
    for k in (1, -1):
        # for i in range(neigs-int((1+k)/2)):
        for i in range(1):
            data = masses[k][i]
            label = "lam={}, g={}".format(lam,g)
            plt.plot(xlist, data, label=label,
                    markersize=markersize, marker=marker)
        plt.gca().set_prop_cycle(None)

argv = sys.argv


if len(argv) < 5:
    print("{} <L> <g4> <Emax> <Lambda>".format(argv[0]))
    sys.exit(-1)

L = float(argv[1])
g4 = float(argv[2])
Emax = float(argv[3])
lam  = float(argv[4])

lammin = 4
ETmin = 10

g4list = np.linspace(0,30,15)
lamlist = np.linspace(lammin, Lambda, 10)
ETlist = np.linspace(ETmin, Emax, 10)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))

plotvsET(L=L, lam=lam, g2=g2, g4=g4, ETlist=ETlist)

title = r"g4={}, L={}, Lambda={}".format(g4, L, Lambda)
fname = r"g4={}, L={}, Lambda={}".format(g4, L, Lambda)
loc = "upper right"

# Vacuum
plt.figure(1)
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_T$")
plt.ylabel(r"$\mathcal{E}_0$")
plt.legend(loc=loc)
plt.savefig("vacvsET"+fname)


# Mass
plt.figure(2)
plt.figure(3, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$m_{\rm ph}$")
plt.legend(loc=loc)
plt.savefig("massvsET"+fname)
