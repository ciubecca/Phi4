import sys
import numpy as np
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
import database
from sys import exit, argv
from time import time
from paramplots import *
from counterterms import *

form  = "pdf"

logct = True
fourfacnorm = False

# Subtract local vacuum counterterm from the plot?
subvac = True
useexactct = True
cubic = True
exactcubic = True

ETmin = 10
nET = 10

klist = (1,-1)
neigs = 4

db = database.Database()

# ETmaxdict = {5:24}
ETmaxdict = {5:24, 6:20, 7:18}


def plotvsET(L, g2, g4, ETlist, ct):

    xlist = ETlist

    spectrum = {k:[] for k in klist}
    masses = {}

    for k in klist:
        for ET in ETlist:

            approxQuery = {"g4":g4, "g2":g2, "L":L, "ET":ET}
            exactQuery = {"k": k, "neigs":neigs, "logct":logct, "impr":False, "fourfacnorm":fourfacnorm}
            boundQuery = {}
            exactQuery["momcut"] = False

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

            if subvac:
                if ct!=None:
                    data -= g4**2*array([ct.ct2(ET, 0) for ET in ETlist])
                else:
                    data -= g4**2*array([ct0ET(ET, 0, 1) for ET in ETlist])

                if cubic:
                    if exactcubic:
                        data -= g4**3*array([ct.ct3(ET) for ET in ETlist])
                    else:
                        data -= g4**3*array([ct0ET3(ET, 1) for ET in ETlist])

            label = r"$L$={}".format(L)
            plt.plot(ETlist, data, label=label, color=color[k])

    # MASS
    plt.figure(2)
    # for k in (-1,1):
    for k in (1, -1):
        # for i in range(neigs-int((1+k)/2)):
        for i in range(1):
            data = masses[k][:,i]
            if k==-1: print("L={}, ET={}, mass={}".format(L, ETlist[0], data[0]))
            label = r"$L$={}, $k$={}".format(L,k)
            plt.plot(xlist, data, label=label, color=color[k])

argv = sys.argv


if len(argv) < 3:
    print("{} <g4> <g2>".format(argv[0]))
    sys.exit(-1)

g4 = float(argv[1])
g2 = float(argv[2])


for i,(L,ETmax) in enumerate(ETmaxdict.items()):

    if useexactct:
        print("Computing matrices for exact counterterms...")
        ct = exactct(L, ETmax)
        print("Done")
    else:
        ct = None

    ETlist = np.linspace(ETmin, ETmax, nET)
    setparams(i)
    plotvsET(L=L, g2=g2, g4=g4, ETlist=ETlist, ct=ct)


title = r"g2={}, g4={}".format(g2, g4)
fname = r"g2={}, g4={}".format(g2, g4)
loc = "upper right"

# Vacuum
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_T$")
plt.ylabel(r"$\mathcal{E}_0/L^2$")
plt.legend(loc=loc)
plt.savefig("plots/vacvsET_{}_{}_{}_{}_{}.{}".format(fname,subvac,useexactct,
    cubic,exactcubic,form))
plt.clf()


# Mass
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_T$")
plt.ylabel(r"$m_{\rm ph}$")
plt.legend(loc=loc)
plt.savefig("plots/massvsET_{}.{}".format(fname,form))
plt.clf()
