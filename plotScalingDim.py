import sys
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
import math
from scipy import pi, log, log10, array, sqrt, stats, exp
from matplotlib import rc
from cycler import cycler
import database
from sys import exit
import numpy as np
from numpy import concatenate as concat
from extrapolate import *

neigs = {}
neigs[1] = 2
neigs[-1] = 3

Llist = [5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]

output = "pdf"
# renlist = ("rentails", "renloc")
renlist = ("rentails", )

marker = 'o'
markersize = 2.5

plt.style.use('ggplot')
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

klist = (1,-1)
color = {}
color[1] = 'b'
color[-1] = 'r'

ymin = {}
ymax = {}
xmax = max(Llist)+0.1
xmin = min(Llist)-0.1

db = database.Database("data/spectra3.db")


def plotvsL(Llist, gmin, gmax):

    Lambda = [np.zeros(len(Llist))]*2
    Spectrum = [{k: np.zeros((neigs[k], len(Llist))) for k in klist}]*2

    for j,g in enumerate(gmin, gmax):
        for i,L in enumerate(Llist):
            e = {}
            e[1] = Extrapolator(db, 1, L, g)
            e[-1] = Extrapolator(db, 1, L, g)

            e.train()

            Lambda[j][i] = E[1][0]/L
            for k in klist:
                if k==1: n=1
                else: n=0
                Spectrum[j][k][:, i] = (E[k][n:neigs[k]+n]-E[1][0])/L


    # Lambda
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    plt.fill_between(Llist, Lambda[0], Lambda[1], marker=marker)

    # Mass
    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    for k in (-1,1):
        for n in range(neigs[k]):
            if n ==0:
                label = "k={}".format(k)
            else:
                label = None
            plt.fill_between(Llist, Spectrum[0][k][n], Spectrum[1][k][n],
                    c = color[k], label=label)

argv = sys.argv

if len(argv) < 3:
    print(argv[0], "<gmin> <gmax>")
    sys.exit(-1)

gmin = float(argv[1])
gmax = float(argv[2])

print("g=", gmin, gmax)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plotvsL(Llist, gmin, gmax)

fname = "{1}".format( output)

# VACUUM ENERGY DENSITY
title = r"$g$=[{}, {}]".format(gmin, gmax)
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$\mathcal{E}_0/L$")
ymargin = 10**(-5)
plt.xlim(xmin, xmax)
plt.legend(loc=2)
plt.savefig("LambdaCritvsL_"+fname, bbox_inches='tight')

# MASS
title = r"$g$=[{}, {}]".format(gmin, gmax)
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$(\mathcal{E}_1-\mathcal{E}_0)/L$")
plt.xlim(xmin, xmax)
plt.legend(loc=1)

plt.savefig("SpecCritvsL_"+fname, bbox_inches='tight')
