import sys
import numpy as np
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from paramplots import *
import database
from sys import exit, argv

form  = "pdf"

logct = True

params = {'legend.fontsize': 8}
plt.rcParams.update(params)


g4max = 15

g4list = np.linspace(1,g4max,15)
print("g4: ", g4list)


lammin = 4
ETmin = 10

nlam = 3
nET = 10

klist = (1,-1)
neigs = 4

color = {1:"b", -1:"r"}


db = database.Database()


def plotvsg(L, lam, g2, g4list, ET):

    xlist = g4list

    spectrum = {k:[] for k in klist}
    masses = {}

    for k in klist:
        for g4 in g4list:

            approxQuery = {"g4":g4, "g2":g2, "L":L, "ET":ET}
            exactQuery = {"k": k, "neigs":neigs, "logct":True, "impr":False}
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


    # MASS
    plt.figure(2)
    for k in (-1,):
    # for k in (1, -1):
        # for i in range(neigs-int((1+k)/2)):
        for i in range(1):
            data = masses[k][:,i]
            label = r"$\Lambda$={}, $k$={}".format(lam,k)
            plt.plot(xlist, data, label=label, color=color[k])

argv = sys.argv


if len(argv) < 4:
    print("{} <L> <ET> <Lambdamax> <g2>".format(argv[0]))
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
Lambdamax  = float(argv[3])
g2 = float(argv[4])

lamlist = np.linspace(lammin, Lambdamax, nlam)
lamlist = np.concatenate([lamlist, [np.inf]])


for idx, lam in enumerate(lamlist):
    setparams(idx)
    plotvsg(L=L, lam=lam, g2=g2, ET=ET, g4list=g4list)

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
plt.title(title)
plt.xlabel(r"$g_4$")
plt.ylabel(r"$m_{\rm ph}$")
plt.legend(loc=loc)
plt.savefig("plots/massvsg_{}.{}".format(fname,form))
plt.clf()

plt.gca().set_prop_cycle(None)
