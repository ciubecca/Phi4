import sys
import matplotlib.pyplot as plt
import matplotlib
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
neigs[1] = 1
neigs[-1] = 2

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


def plotvsL(Llist, g):

    Spectrum = {k: np.zeros((neigs[k], len(Llist))) for k in klist}
    Lambda = np.zeros(len(Llist))

    LambdaErr = np.zeros((2, len(Llist)))

    for i,L in enumerate(Llist):
#           print("L:{}".format(L))
        e = {}
        e[1] = Extrapolator(db, 1, L, g)
        e[-1] = Extrapolator(db, -1, L, g)

        e[1].train(neigs=neigs[1]+1)
        e[-1].train(neigs=neigs[-1])

        Lambda[i] = e[1].asymValue()[0]/L

        print(e[1].asymErr().shape)
        # LambdaErr[:,i] =

        for k in klist:
            if k==1: n=1
            else: n=0
            Spectrum[k][:, i] = (e[k].asymValue()[n:neigs[k]+n]-
                    e[1].asymValue()[0])*L/(2*pi)

    # Lambda
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    plt.fill_between(Llist, Lambda[0], Lambda[1])

    # Mass
    fig = plt.figure(2)
    ax = fig.add_subplot(111)

    trans = matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData)

    plt.axhline(y=1/8, color='k', linestyle='--', dashes=[10,10], linewidth=0.9)
    ax.annotate(r"$\Delta_\sigma$", xy=(1.01,1/8), xycoords=trans, clip_on=False,
            va='center')

    plt.axhline(y=2+1/8, color='k', linestyle='--', dashes=[10,10], linewidth=0.9)
    ax.annotate(r"$\Delta_{\partial^2 \sigma}$",
            xy=(1.01,2+1/8), xycoords=trans, clip_on=False, va='center')


    plt.axhline(y=1, color='k', linestyle='--', dashes=[10,10], linewidth=0.9)
    ax.annotate(r"$\Delta_\epsilon$", xy=(1.01,1), xycoords=trans, clip_on=False,
            va='center')

    for k in (-1,1):
        for n in range(neigs[k]):
            if n ==0:
                label = "k={}".format(k)
            else:
                label = None
            plt.fill_between(Llist, Spectrum[0][k][n], Spectrum[1][k][n],
                    color=color[k], label=label)

argv = sys.argv

if len(argv) < 2:
    print(argv[0], "<g>")
    sys.exit(-1)

g = float(argv[1])

print("g=", g)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plotvsL(Llist, g)

fname = "{}".format(output)

# VACUUM ENERGY DENSITY
title = r"$g$={}".format(g)
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$\mathcal{E}_0/L$")
ymargin = 10**(-5)
plt.xlim(xmin, xmax)
plt.legend(loc=2)
plt.savefig("LambdaCritvsL."+fname, bbox_inches='tight')

# MASS
title = r"$g$={}".format(g)
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$(\mathcal{E}_I-\mathcal{E}_0) L/(2 \pi)$")
plt.xlim(xmin, xmax)
plt.legend(loc=1)

plt.savefig("SpecCritvsL."+fname, bbox_inches='tight')
