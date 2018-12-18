import sys
import numpy as np
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
import database
from paramplots import *
from sys import exit, argv

form  = "pdf"

g4max = 30

Llist = [5, 6, 7]
ETlist = [24, 20, 18]

Llist = [5, 7]
ETlist = [24, 19]

g4list = np.linspace(1, g4max, 30)
print("g4: ", g4list)


klist = (1,-1)
neigs = 4

db = database.Database()


def plotvsg(L, lam, g2, g4list, ET):

    xlist = g4list

    spectrum = {k:[] for k in klist}
    masses = {}

    for k in klist:
        for g4 in g4list:

            approxQuery = {"g4":g4, "g2":g2, "L":L, "ET":ET, "Lambda":lam}
            exactQuery = {"k": k, "neigs":neigs, "logct":True}
            boundQuery = {}


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
            label = r"$E_T$={}, $L$={}".format(ET,L)
            plt.plot(xlist, data, label=label, color=color[k])


    # MASS
    plt.figure(2)
    for k in (-1,):
    # for k in (1, -1):
        # for i in range(neigs-int((1+k)/2)):
        for i in range(1):
            data = masses[k][:,i]
            label = r"$E_T$={}, $L$={}".format(ET,L)
            plt.plot(xlist, data, label=label, color=color[k])

argv = sys.argv


if len(argv) < 3:
    print("{}  <Lambda> <g2>".format(argv[0]))
    sys.exit(-1)

Lambda  = float(argv[1])
g2 = float(argv[2])


for idx, (L, ET) in enumerate(zip(Llist,ETlist)):
    setparams(idx)
    plotvsg(L=L, lam=Lambda, g2=g2, ET=ET, g4list=g4list)

title = r"$g_2$={}, $\Lambda$={}".format(g2, Lambda)
fname = r"g2={}, Lambda={}".format(g2, Lambda)
loc = "upper right"

# Vacuum
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$g_4$")
plt.ylabel(r"$\mathcal{E}_0/L$")
plt.legend(loc=loc)
plt.savefig("plots/vacvsgL_{}.{}".format(fname,form))
plt.clf()


# Mass
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$g_4$")
plt.ylabel(r"$m_{\rm ph}$")
plt.legend(loc=loc)
plt.savefig("plots/massvsgL_{}.{}".format(fname,form))
plt.clf()
