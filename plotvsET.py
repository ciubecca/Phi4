import sys
import numpy as np
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
import database
from sys import exit, argv
from time import time
from paramplots import *

form  = "pdf"

ct = True


lammin = 4
ETmin = 10
nlam = 3
nET = 10

klist = (1,-1)
neigs = 4

db = database.Database()


def plotvsET(L, lam, g2, g4, ETlist):

    xlist = ETlist

    spectrum = {k:[] for k in klist}
    masses = {}

    for k in klist:
        for ET in ETlist:

            approxQuery = {"g4":g4, "g2":g2, "L":L, "ET":ET}
            exactQuery = {"k": k, "neigs":neigs, "logct":ct, "impr":False}
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
            data = spectrum[k][:,i]/L**2
            label = r"$\Lambda$={}".format(lam,g4)
            plt.plot(ETlist, data, label=label, color=color[k])

    # MASS
    plt.figure(2)
    # for k in (-1,1):
    for k in (1, -1):
        # for i in range(neigs-int((1+k)/2)):
        for i in range(1):
            data = masses[k][:,i]
            label = r"$\Lambda$={}, $k$={}".format(lam,k)
            plt.plot(xlist, data, label=label, color=color[k])

argv = sys.argv


if len(argv) < 6:
    print("{} <L> <ETmax> <Lambdamax> <g4> <g2>".format(argv[0]))
    sys.exit(-1)

L = float(argv[1])
ETmax = float(argv[2])
Lambdamax  = float(argv[3])
g4  = float(argv[4])
g2  = float(argv[5])

lamlist = np.linspace(lammin, Lambdamax, nlam)
lamlist = np.concatenate([lamlist, [np.inf]])

ETlist = np.linspace(ETmin, ETmax, nET)


for i,lam in enumerate(lamlist):
    setparams(i)
    plotvsET(L=L, lam=lam, g2=g2, g4=g4, ETlist=ETlist)

title = r"g2={}, g4={}, L={}".format(g2, g4, L)
fname = r"g2={}, g4={}, L={}".format(g2, g4, L)
loc = "upper right"

# Vacuum
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_T$")
plt.ylabel(r"$\mathcal{E}_0/L$")
plt.legend(loc=loc)
plt.savefig("plots/vacvsET_{}.{}".format(fname,form))
plt.clf()


# Mass
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_T$")
plt.ylabel(r"$m_{\rm ph}$")
plt.legend(loc=loc)
plt.savefig("plots/massvsET_{}.{}".format(fname,form))
plt.clf()
