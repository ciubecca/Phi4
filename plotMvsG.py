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

nparams = 3

plt.style.use('ggplot')
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

klist = (1,-1)

glist = scipy.linspace(0.2, 2.4, 12)

xmax = max(glist)+0.01
xmin = min(glist)-0.01

db = database.Database("data/spectra3.db")


def plotvsG(Llist, axes):

    MassInf = np.zeros(len(glist))
    MassErr = np.zeros(len(glist))

    for i,g in enumerate(glist):
        a = ExtrvsL(db, g)
        a.train(nparam=3)
        MassInf[i] = a.asymValue()[-1]
        MassErr[i] = a.asymErr()[-1]

    # Plot extrapolated data
    axes.errorbar(glist, MassInf, MassErr, label=r"$E_T=\infty$")
    print(glist)
    print("{"+",".join("{:f}".format(x) for x in MassInf)+"}")
    print("{"+",".join("{:f}".format(x) for x in MassErr)+"}")

argv = sys.argv

if len(argv) < 1:
    print(argv[0], "")
    sys.exit(-1)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

f, axes = plt.subplots(1, 1)
f.subplots_adjust(hspace=0, wspace=0, top=0.94, right=0.95, left=0)

plotvsG(glist, axes)

axes.set_ylabel(r"$m_{ph}$")
axes.set_ylim(-0.01, 1)
axes.set_xlabel(r"$g$")

fname = ".{0}".format(output)

# MASS
plt.figure(1, figsize=(4, 2.5), dpi=300)
s = "MvsGraw"
if nparams==2:
    s += "_p=2"
plt.savefig(s+fname, bbox_inches='tight')
