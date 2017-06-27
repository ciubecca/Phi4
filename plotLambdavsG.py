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
from matplotlib.ticker import MaxNLocator
from numpy import concatenate as concat
import itertools
from extrapolate import *
from sklearn.linear_model import Ridge, LinearRegression

output = "pdf"
# renlist = ("rentails", "renloc")
renlist = ("rentails", )

labelsize = 20
nparams = 3
L = 10

plt.style.use('ggplot')
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

klist = (1,-1)

glist = scipy.linspace(0.2, 2.6, 13)
glistTot = scipy.linspace(0.2, 3, 15)

xmax = max(glistTot)+0.01
xmin = 0

db = database.Database("data/spectra3.db")


def plotvsG(glist, axes):

    LambdaFiniteL = np.zeros(len(glistTot))
    LambdaErrFiniteL = np.zeros((2, len(glistTot)))
    LambdaInf = np.zeros(len(glist))
    LambdaErr = np.zeros(len(glist))

    for i,g in enumerate(glist):
        a = ExtrvsL(db, g)
        a.train(nparam=3)
        LambdaInf[i] = a.asymValue(1)
        LambdaErr[i] = a.asymErr(1)

    for i,g in enumerate(glistTot):
        b = Extrapolator(db, L=L, g=g)
        b.train(neigs=1)
        LambdaFiniteL[i] = b.asymValue(1)[0]/L
        LambdaErrFiniteL[:,i] = b.asymErr(1)[0]/L

    # Plot extrapolated data
    axes[0].errorbar(glistTot, LambdaFiniteL, LambdaErrFiniteL,
            label=r"$L = {}$".format(L), color="green", ls='none', capthick=1.5)
    axes[0].errorbar(glist, LambdaInf, LambdaErr, label=r"$L=\infty$", ls="--",
            capthick=1.5)

    # Plot just error bars, divided by g^2
    axes[1].errorbar(glist, [0]*len(glist),
            LambdaErr/glist**2, ls='none')

    print("glist:", glist)
    print("LambdaInf:")
    print("{"+",".join("{:f}".format(x) for x in LambdaInf)+"}")
    print("LambdaErr:")
    print("{"+",".join("{:f}".format(x) for x in LambdaErr)+"}")

    print("LambdaFiniteL:")
    print("{"+",".join("{:f}".format(x) for x in LambdaFiniteL)+"}")
    print("LambdaErrFiniteL:")
    print("{"+",".join("{"+"{0:f}, {1:f}".format(x[0],x[1])+"}"
        for x in LambdaErrFiniteL.transpose())+"}")

argv = sys.argv

if len(argv) < 1:
    print(argv[0], "")
    sys.exit(-1)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

f, axes = plt.subplots(2, 1, sharex='col', gridspec_kw = {'height_ratios':[3,1]})
f.subplots_adjust(hspace=0, wspace=0, top=0.94, right=0.95, left=0)

plotvsG(glist, axes)

axes[0].set_ylabel(r"$\mathcal{E}_0/L$", fontsize=labelsize)
# axes[0].set_ylim(-0.01, 1.01)
axes[1].set_ylabel(r"errors$/g^2$", fontsize=labelsize)
axes[1].set_xlabel(r"$g$", fontsize=labelsize)
axes[1].set_xlim(xmin, xmax)
axes[0].legend(fontsize=12)
# axes[1].set_ylim(-0.03, 0.03)
# axes[1].yaxis.set_ticks(np.arange(-0.001, 0.001, 5))
axes[1].locator_params(nbins=8, axis='y')
axes[1].yaxis.set_major_locator(MaxNLocator(prune='upper', nbins=6))


fname = ".{0}".format(output)

# LAMBDA
plt.figure(1, figsize=(4, 2.5), dpi=300)
s = "LambdavsG"
if nparams==2:
    s += "_p=2"
plt.savefig(s+fname, bbox_inches='tight')
