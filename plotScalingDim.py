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

Llist = np.array([5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10])

output = "pdf"

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
    SpectrumErr = {k: np.zeros((neigs[k], 2, len(Llist))) for k in klist}

    Lambda = np.zeros(len(Llist))
    LambdaErr = np.zeros((2, len(Llist)))

    for i,L in enumerate(Llist):
        e = {}
        e[1] = Extrapolator(db, 1, L, g)
        e[-1] = Extrapolator(db, -1, L, g)

        e[1].train(neigs=neigs[1]+1)
        e[-1].train(neigs=neigs[-1])

        Lambda[i] = e[1].asymValue()[0]/L
        LambdaErr[:,i] = e[1].asymErr()[0]/L

        for k in klist:
            if k==1: n=1
            else: n=0
            Spectrum[k][:, i] = (e[k].asymValue()[n:neigs[k]+n]-
                    e[1].asymValue()[0])*L/(2*pi)

            # Error bars of the eigenvalues
            en = e[k].asymErr()[n:neigs[k]+n]
            # Error bars of the vacuum to subtract
            e0 = np.repeat(e[1].asymErr()[[0],::-1], neigs[k], axis=0)

            # XXX
            SpectrumErr[k][:, :, i] = np.amax(np.array([e0, en]), axis=0)*L/(2*pi)
            SpectrumErr[k][:, :, i] = (e0+en)*L/(2*pi)

    # Lambda
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    plt.errorbar(Llist, Lambda, LambdaErr)


    X = np.array([1/Llist**2]).transpose()
    Y = Lambda
    model = LinearRegression().fit(X, Y, sample_weight=1/np.amax(LambdaErr,axis=0))
    xlist = np.linspace(min(Llist), max(Llist), 100)
    plt.plot(xlist, model.predict(np.array([1/xlist**2]).transpose()))

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
            plt.errorbar(Llist, Spectrum[k][n], SpectrumErr[k][n],
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
