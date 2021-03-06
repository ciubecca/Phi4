##########################################
# Plots dependence of spectrum on g for a given L, at Lambda=infinity (no Lambda cutoff)
# Also, add Chang duality predictions to the plot
##########################################

import sys
import numpy as np
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from paramplots import *
import database
from sys import exit, argv
from chang import *

form  = "pdf"

logct = True
fourfacnorm = False

params = {'legend.fontsize': 8}
plt.rcParams.update(params)


g4list = np.linspace(0.2,6,30)
print("g4: ", g4list)

# Index of the value of g in glist which is closest to gstar
istarapprox = np.argmin(abs(g4list-gstar))
print(istarapprox)


ETmin = 10
nET = 10

klist = (1,-1)
neigs = 4

color = {1:"b", -1:"r"}


db = database.Database()


def plotvsg(L, g2, g4list, ET):

    lam = np.inf
    xlist = g4list

    spectrum = {k:[] for k in klist}
    masses = {}

    massesDual = {}
    gdual = {}


    for k in klist:
        for g4 in g4list:

            approxQuery = {"g4":g4, "g2":g2, "L":L, "ET":ET}
            exactQuery = {"k": k, "neigs":neigs, "logct":logct, "fourfacnorm":fourfacnorm, "impr":False}
            boundQuery = {}

            if lam==np.inf:
                exactQuery["momcut"] = False
            else:
                exactQuery["momcut"] = True
                exactQuery["Lambda"] = lam

            try:
                spectrum[k].append(db.getObjList('spec', exactQuery, approxQuery, orderBy="date")[0])

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
        masses[k] = (spectrum[k][:,imin:].transpose()-spectrum[1][:,0]).transpose()


    # Compute Chang duality predictions

    # Impose lower bound on g so that the dual g is not too large
    dualidxlist = np.array([i for i,g in enumerate(g4list) if g>1.0 and g<gstar])
    gFirstBranch = g4list[dualidxlist]
    gDual = array([xmintomax2(g) for g in gFirstBranch])
    # Scale factor to go from first to second branch
    factordual = array([factorToSecondBranch2(g) for g in gFirstBranch])

    # Dual lowest mass gap
    # massesDual = masses[-1][dualidxlist,0]
    massesDual = masses[-1][dualidxlist,0]*factordual

    print("gFirstBranch:", gFirstBranch)
    print("gDual:", gDual)
    print("factorDual:", factordual)
    print("massesDual:", massesDual)


    # MASS
    plt.figure(2)
    for k in (-1,):
        for i in range(1):
            data = masses[k][:,i]

            label = "Chang, $E_T={:.4f}$".format(ET)
            plt.plot(gDual, massesDual, label=label, color='b', linewidth=1)

            label = r"$E_T={:.4f}$".format(ET)
            plt.plot(xlist, data, label=label, color='r', linewidth=1)

            exstar = data[istarapprox]
            x = np.linspace(gstar, 6, 100)
            plt.plot(x, exstar*x/gstar, linestyle='dashed', c='g',
                    label='analytic')
            # plt.plot(x, exstar*(2*x-gstar)/gstar, linestyle='dashed', c='g',
                    # label='not analytic')

argv = sys.argv


if len(argv) < 3:
    print("{} <L> <ETmax> <g2>".format(argv[0]))
    sys.exit(-1)

L = float(argv[1])
ETmax = float(argv[2])
g2 = float(argv[3])


ETlist = np.linspace(ETmin, ETmax, nET)
print("ETlist: {}".format(ETlist))
ETlistred = ETlist[1::4]
print("ETlist reduced: {}".format(ETlistred))

for i,ET in enumerate(ETlistred):
    setparams(i)
    plotvsg(L=L, g2=g2, ET=ET, g4list=g4list)

title = r"g2={}, ET={}, L={}".format(g2, ET, L)
fname = r"g2={}, ET={}, L={}".format(g2, ET, L)
loc = "upper right"

# Mass
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')

plt.axvline(gstar, c='k', linestyle='dashed', linewidth=1)

plt.xlim(0,6)
plt.ylim(0,5)

plt.title(title)
plt.xlabel(r"$g_4$")
plt.ylabel(r"$m_{\rm ph}$")
plt.legend(loc=loc)
plt.savefig("plots/massvsgChang_{}.{}".format(fname,form))
plt.clf()

plt.gca().set_prop_cycle(None)
