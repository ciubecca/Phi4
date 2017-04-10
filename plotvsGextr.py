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

output = "png"
# renlist = ("rentails", "renloc")
renlist = ("rentails", )


plt.style.use('ggplot')
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

klist = (1,-1)

glist = scipy.linspace(0.2, 3.0, 15)

xmax = max(glist)+0.01
xmin = min(glist)-0.01

db = database.Database("data/spectra3.db")


def plotvsG(Llist):

    MassInf = np.zeros(len(glist))
    MassErr = np.zeros((2, len(glist)))

    for i,g in enumerate(glist):
        a = ExtrvsL(db, g)
        a.train()
        MassInf[i] = a.asymValue()[-1]
        MassErr[:,i] = a.asymErr()[-1]


    print(MassInf)
    print(MassErr)

    # Plot extrapolated data
    plt.errorbar(glist, MassInf, MassErr, label=r"$E_T=\infty$")


argv = sys.argv

if len(argv) < 1:
    print(argv[0], "")
    sys.exit(-1)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plotvsG(glist)

fname = ".{0}".format(output)

# MASS
title = r"".format()
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$g$")
plt.ylabel(r"$m_{ph}$")
plt.xlim(xmin, xmax)
plt.ylim(-0.01, 1)
plt.legend(loc=1)
plt.savefig("MvsG_"+fname, bbox_inches='tight')
