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

form  = "png"

logct = True
fourfacnorm = False

params = {'legend.fontsize': 8}
plt.rcParams.update(params)


g4list = np.linspace(0.2,6,30)
print("g4: ", g4list)


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


    # VACUUM
    plt.figure(1)
    for k in (1,):
        # for i in range(neigs):
        for i in range(1):
            data = spectrum[k][:,i]/L
            label = r"$\Lambda$={}".format(lam,g4)
            plt.plot(xlist, data, label=label, color=color[k])

    # Compute Change duality predictions

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
    # for k in (1, -1):
        # for i in range(neigs-int((1+k)/2)):
        for i in range(1):
            data = masses[k][:,i]
            label = r"$\Lambda$={}, $k$={}".format(lam,k)
            plt.plot(xlist, data, label=label, color=color[k], markersize=3)

            label = "Chang"
            plt.plot(gDual, massesDual, label=label, linestyle='dotted', color='b',
                    markersize=3)

argv = sys.argv


if len(argv) < 3:
    print("{} <L> <ET> <g2>".format(argv[0]))
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
g2 = float(argv[3])

setparams(0)
plotvsg(L=L, g2=g2, ET=ET, g4list=g4list)

title = r"g2={}, ET={}, L={}".format(g2, ET, L)
fname = r"g2={}, ET={}, L={}".format(g2, ET, L)
loc = "upper right"

# Vacuum
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$g_4$")
plt.ylabel(r"$\mathcal{E}_0/L$")
plt.legend(loc=loc)
plt.savefig("plots/vacvsg_{}.{}".format(fname,form))
plt.clf()


# Mass
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')

plt.axvline(gstar)

plt.title(title)
plt.xlabel(r"$g_4$")
plt.ylabel(r"$m_{\rm ph}$")
plt.legend(loc=loc)
plt.savefig("plots/massvsg_{}.{}".format(fname,form))
plt.clf()

plt.gca().set_prop_cycle(None)
