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

gdratio = 1.618
figsize = 5
labelsize=20

neigs = 2

Llist = np.array([5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10])
Llist = np.array([6,6.5,7,7.5,8,8.5,9,9.5,10])

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


def plotvsL(Llist, gpair):

    Spectrum = {g: {k: np.zeros((neigs, len(Llist))) for k in klist} for g in gpair}
    SpectrumMax = {g: {k: np.zeros((neigs, len(Llist))) for k in klist} for g in gpair}
    SpectrumMin = {g: {k: np.zeros((neigs, len(Llist))) for k in klist} for g in gpair}
    SpectrumErr = {g: {k: np.zeros((neigs, 2, len(Llist))) for k in klist} for g in gpair}

    Lambda = {g: np.zeros(len(Llist)) for g in gpair}
    LambdaErr = {g: np.zeros((2, len(Llist))) for g in gpair}

    for g in gpair:
        for i,L in enumerate(Llist):
            e = Extrapolator(db, L, g)
            e.train(neigs=neigs)

            for k in klist:
                Spectrum[g][k][:, i] = e.asymValue(k)*L/(2*pi)
                Spectrum[g][1][0, i] = e.asymValue(1)[0]

                # Error bars of the eigenvalues
                SpectrumErr[g][k][:, :, i] = e.asymErr(k)*L/(2*pi)
                SpectrumErr[g][1][:, 0, i] = e.asymErr(1)[:,0]


    for k in (-1,1):
        SpectrumMax[k] = np.amax(np.array(
            [Spectrum[gpair[0]][k]+SpectrumErr[gpair[0]][k][1],
            Spectrum[gpair[1]][k]+SpectrumErr[gpair[1]][k][1]]),
            axis=0)

        SpectrumMin[k] = np.amin(np.array(
            [Spectrum[gpair[0]][k]-SpectrumErr[gpair[0]][k][0],
            Spectrum[gpair[1]][k]-SpectrumErr[gpair[1]][k][0]]),
            axis=0)


    # Lambda
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    for g in gpair:
        err = SpectrumErr[g][1][:,0]/Llist
        Y = Spectrum[g][1][0]/Llist
        plt.errorbar(Llist, Y, err)

    # X = np.array([1/Llist**2]).transpose()
    # model = LinearRegression().fit(X, Y, sample_weight=1/np.amax(err,axis=0))
    # xlist = np.linspace(min(Llist), max(Llist), 100)
    # plt.plot(xlist, model.predict(np.array([1/xlist**2]).transpose()))

    # Mass
    fig = plt.figure(2)
    ax = fig.add_subplot(111)

    trans = matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData)

    plt.axhline(y=1/8, color='k', linestyle='--', dashes=[10,10], linewidth=0.9)
    ax.annotate(r"$\Delta_\sigma$", xy=(1.01,1/8), xycoords=trans, clip_on=False,
            va='center', fontsize=15)

    plt.axhline(y=2+1/8, color='k', linestyle='--', dashes=[10,10], linewidth=0.9)
    ax.annotate(r"$\Delta_{\partial^2 \sigma}$", fontsize=15,
            xy=(1.01,2+1/8), xycoords=trans, clip_on=False, va='center')


    plt.axhline(y=1, color='k', linestyle='--', dashes=[10,10], linewidth=0.9)
    ax.annotate(r"$\Delta_\epsilon$", xy=(1.01,1), xycoords=trans, clip_on=False,
            va='center', fontsize=15)

    # for g in gpair:
    for k in (-1,1):
        for n in range(neigs):
            if n ==0 and k==1:
                continue
            elif (n==0 and k==-1) or (n==1 and k==1):
                label = "k={}".format(k)
            else:
                label = None
            # plt.errorbar(Llist, Spectrum[g][k][n], SpectrumErr[g][k][n],
                    # color=color[k], label=label)

            plt.fill_between(Llist, SpectrumMax[k][n], SpectrumMin[k][n],
                    color="r")

argv = sys.argv

if len(argv) < 2:
    print(argv[0], "<gmin> <gmax>")
    sys.exit(-1)

gpair = []
gpair.append(float(argv[1]))
gpair.append(float(argv[2]))

print("g={}".format(gpair))

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plotvsL(Llist, gpair)

fname = "{}".format(output)

# VACUUM ENERGY DENSITY
title = r"$g \in $ [{},{}]".format(gpair[0],gpair[1])
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$\mathcal{E}_0/L$")
ymargin = 10**(-5)
plt.xlim(xmin, xmax)
plt.legend(loc=2)
plt.savefig("LambdaCritvsL."+fname, bbox_inches='tight')

# MASS
fig = plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
fig.set_figheight(figsize)
fig.set_figwidth(figsize*gdratio)
plt.title(title)
plt.xlabel(r"$L$", fontsize=labelsize)
plt.ylabel(r"$(\mathcal{E}_I-\mathcal{E}_0) L/(2 \pi)$", fontsize=labelsize)
plt.xlim(xmin, xmax)
# plt.legend(loc=1, fontsize=12)

plt.savefig("SpecCritvsL."+fname, bbox_inches='tight')
