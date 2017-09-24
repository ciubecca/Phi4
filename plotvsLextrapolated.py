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

nparam = 3
gdratio = 1.618
figsize = 5
labelsize = 15

plt.rc('axes.formatter', useoffset=False)

Llist = {}
Llist["rentails"] = [5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]
Llist["renloc"] = [5,5.5,6,6.5,7,7.5,8,8.5,9,10]

label = {"rentails":"NLO", "renloc":"local LO"}

output = "pdf"
renlist = ("rentails", "renloc")
renlist = ("rentails",)

marker = '*'
markersize = 2.5

plt.style.use('ggplot')
plt.rc('axes', prop_cycle=(
    cycler('marker', ['*', '', 'o', '^'])
    +cycler('color', ['r', 'g', 'b', 'y'])
    ))


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

klist = (1,-1)

ymin = {}
ymax = {}
xmax = max(Llist["rentails"])+0.1
xmin = min(Llist["rentails"])-0.1

db = database.Database("data/spectra3.db")


def plotvsL(Llist):

    Lambda = {ren: np.zeros(len(Llist[ren])) for ren in renlist}
    Mass = {ren: np.zeros(len(Llist[ren])) for ren in renlist}

    for ren in renlist:
        for i,L in enumerate(Llist[ren]):
            E0 = db.getEigs(1, ren, g, L, ETmax[ren][L])[0]
            E1 = db.getEigs(-1, ren, g, L, ETmax[ren][L])[0]
            Lambda[ren][i] = E0/L
            Mass[ren][i] = (E1-E0)

    a = ExtrvsL(db, g)
    a.train(nparam)

    ymax[1] = max(max(max(Lambda[ren]) for ren in renlist),
            max(a.LambdaInf+a.LambdaErr[1]))
    ymin[1] = min(min(min(Lambda[ren]) for ren in renlist),
            min(a.LambdaInf-a.LambdaErr[0]))
    ymax[-1] = max(max(max(Mass[ren]) for ren in renlist),
            max(a.MassInf+a.MassErr[1]))
    ymin[-1] = min(min(min(Mass[ren]) for ren in renlist),
            min(a.MassInf-a.MassErr[0]))

    # Lambda
    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    # Plot extrapolated data

    ax.errorbar(Llist["rentails"], a.LambdaInf, a.LambdaErr,
            marker=marker, label=r"$E_T=\infty$")

    # Plot fitted function
    xdata = scipy.linspace(xmin, xmax, 100)
    ax.plot(xdata, a.predict(1, xdata), label="fit")

    # ax.text(0.7, 0.1, "f(L) = "+fvacStr, horizontalalignment='center',
            # verticalalignment='center', fontsize=12, transform=ax.transAxes)

    # for i, m in enumerate(a.msg[1]):
        # ax.text(0.8, 0.1-i*0.05, a.msg[1][i], horizontalalignment='center',
            # verticalalignment='center', fontsize=13, transform=ax.transAxes)

    # Plot non-extrapolated data
    for ren in renlist:
        plt.plot(Llist[ren], Lambda[ren], label=label[ren])


    # Mass
    fig = plt.figure(2)
    ax = fig.add_subplot(111)

    # Plot extrapolated data
    ax.errorbar(LList, a.MassInf, a.MassErr, marker=marker, label=r"$E_T=\infty$")

    # Plot fits
    xdata = scipy.linspace(xmin, xmax, 100)
    ax.plot(xdata, a.predict(-1, xdata), label="fit")

    # Print fitted parameters
    # for i, m in enumerate(a.msg[-1]):
        # ax.text(0.8, 0.85-i*0.05, a.msg[-1][i], horizontalalignment='center',
            # verticalalignment='center', fontsize=13, transform=ax.transAxes)

    # Plot non-extrapolated data
    for ren in renlist:
        plt.plot(Llist[ren], Mass[ren], label=label[ren])


argv = sys.argv

if len(argv) < 2:
    print(argv[0], "<g>")
    sys.exit(-1)

g = float(argv[1])

print("g=", g)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plotvsL(Llist)

fname = "g={0:.1f}.{1}".format(g, output)

plt.figure(figsize=(figsize*gdratio,figsize))


# VACUUM ENERGY DENSITY
title = r"$g$={:.1f}, \quad f(x)={}".format(g, fvacStr)
title = r"$g$={:.1f}".format(g)
plt.figure(1,
        # figsize=(figsize*gdratio,figsize),
        dpi=300, facecolor='w', edgecolor='w')
plt.title(title, fontsize=15)
plt.xlabel(r"$L$", fontsize=labelsize)
plt.ylabel(r"$\mathcal{E}_0/L$", fontsize=labelsize)
ymargin = 10**(-5)
plt.xlim(xmin, xmax)
plt.ylim(ymin[1]-ymargin, ymax[1]+ymargin)
plt.legend(loc=2, fontsize=10)
plt.savefig("extrLambdavsL_"+fname, bbox_inches='tight')

# MASS
if nparam == 2:
    f = fmassStr2
else:
    f = fmassStr
title = r"$g$={:.1f}, \quad f(x)={}".format(g, f)
title = r"$g$={:.1f}".format(g)
fig = plt.figure(2,
        # figsize=(figsize*gdratio,figsize),
        dpi=300, facecolor='w', edgecolor='w')
# fig.set_figheight(figsize)
# fig.set_figwidth(figsize*gdratio)
plt.title(title, fontsize=15)
plt.xlabel(r"$L$", fontsize=labelsize)
plt.ylabel(r"$\mathcal{E}_1-\mathcal{E}_0$", fontsize=labelsize)
plt.xlim(xmin, xmax)
plt.ylim(ymin[-1]-ymargin, ymax[-1]+ymargin)
plt.legend(loc=1, fontsize=10)

s = "extrMassvsL_"
if nparam==2:
    s += "p=2_"
plt.savefig(s+fname, bbox_inches='tight')
# plt.savefig(s+fname)
