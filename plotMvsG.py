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

gdratio = 1.618
figsize = 5
labelsize = 20

nparams = 3

plt.style.use('ggplot')
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))

# Imported from Mathematica
gc = 1/0.363112
def fit(g):
    return ((1 - 0.363112*g)*(1 + 3.44561*g + 1.4152*g**2))/\
            ((1 + 0.325964*g)*(1 + 2.75653*g))

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

klist = (1,-1)

glist = scipy.linspace(0.2, 2.6, 13)
glistTot = scipy.linspace(0.2, 3, 15)

# glist = scipy.linspace(0.2, 0.4, 2)
# glistTot = scipy.linspace(0.2, 0.6, 3)

xmax = max(glistTot)+0.03
xmin = 0

db = database.Database("data/spectra3.db")


def plotvsG(Llist, axes):

    MassFiniteL = np.zeros(len(glistTot))
    MassErrFiniteL = np.zeros((2, len(glistTot)))
    MassInf = np.zeros(len(glist))
    MassErr = np.zeros(len(glist))

    for i,g in enumerate(glist):
        a = ExtrvsL(db, g)
        a.train(nparam=3)
        MassInf[i] = a.asymValue(-1)
        MassErr[i] = a.asymErr(-1)

    for i,g in enumerate(glistTot):
        b = Extrapolator(db, L=10, g=g)
        b.train(neigs=1)
        MassFiniteL[i] = b.asymValue(-1)[0]
        MassErrFiniteL[:,i] = b.asymErr(-1)[0]

    # Plot extrapolated data
    axes[0].errorbar(glist, MassInf, MassErr, label=r"$L=\infty$",
            capthick=1.5, ls='none')
    axes[0].errorbar(glistTot, MassFiniteL, MassErrFiniteL,
            label=r"$L = 10$", color="green", ls='none', capthick=1.5)

# Plot fitted function imported from Mathematica
    xlist = scipy.linspace(0, gc, 100)
    axes[0].plot(xlist, fit(xlist), color='blue', label="fit")


# Plot residuals, divided by g^2
    axes[1].errorbar(glist, (MassInf-fit(glist))/glist**2,
            MassErr/glist**2, label=r"$L=\infty$", ls='none')
    # axes[1].errorbar(glist, MassFiniteL[:glist.size]-fit(glist),
            # MassErrFiniteL[:, :glist.size], label=r"$L = 10$", color="green",
            # ls='none')


    print("glist:", glist)
    print("MassInf:")
    print("{"+",".join("{:f}".format(x) for x in MassInf)+"}")
    print("MassErr:")
    print("{"+",".join("{:f}".format(x) for x in MassErr)+"}")

    print("MassFiniteL:")
    print("{"+",".join("{:f}".format(x) for x in MassFiniteL)+"}")
    print("MassErrFiniteL:")
    print("{"+",".join("{"+"{0:f}, {1:f}".format(x[0],x[1])+"}"
        for x in MassErrFiniteL.transpose())+"}")

argv = sys.argv

if len(argv) < 1:
    print(argv[0], "")
    sys.exit(-1)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

f, axes = plt.subplots(2, 1, sharex='col',
        # figsize=(figsize*gdratio,figsize),
        gridspec_kw = {'height_ratios':[3,1]}
        )
f.subplots_adjust(hspace=0, wspace=0, top=0.94, right=0.95, left=0)

axes[0].set_ylabel(r"$m_{\rm ph}$", fontsize=labelsize)
axes[0].set_ylim(-0.01, 1.01)
axes[1].set_ylabel(r"residuals$/g^2$", fontsize=labelsize)
axes[1].set_xlabel(r"$g$", fontsize=labelsize)
axes[1].set_xlim(xmin, xmax)
# axes[1].set_ylim(-0.03, 0.03)
# axes[1].yaxis.set_ticks(np.arange(-0.001, 0.001, 5))
axes[1].locator_params(nbins=8, axis='y')
axes[1].yaxis.set_major_locator(MaxNLocator(prune='upper', nbins=6))

plotvsG(glist, axes)
axes[0].legend(loc=1, fontsize=12)

fname = ".{0}".format(output)

# MASS
plt.figure(1, figsize=(4, 2.5), dpi=300)
s = "MvsG"
if nparams==2:
    s += "_p=2"
plt.savefig(s+fname, bbox_inches='tight')
