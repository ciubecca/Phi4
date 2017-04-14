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

glist = scipy.linspace(0.2, 2.6, 13)

xmax = max(glist)+0.01
xmin = min(glist)-0.01

db = database.Database("data/spectra3.db")

gc = 0

def mfun(g, gc, a):
    return a*(gc-g)

def plotvsG(Llist, axes):

    global gc

    MassInf = np.zeros(len(glist))
    MassErr = np.zeros(len(glist))
    MassInf2 = np.zeros(len(glist))
    MassErr2 = np.zeros(len(glist))

    for i,g in enumerate(glist):
        a = ExtrvsL(db, g)
        a.train(nparam=3)
        MassInf[i] = a.asymValue()[-1]
        MassErr[i] = a.asymErr()[-1]

        a = ExtrvsL(db, g)
        a.train(nparam=2)
        MassInf2[i] = a.asymValue()[-1]
        MassErr2[i] = a.asymErr()[-1]


    # Plot extrapolated data
    axes[0].errorbar(glist, MassInf, MassErr, label=r"$E_T=\infty$")

    mask = np.logical_and(1.3<glist,  glist<2.5)
    xfit = glist[mask]
    yfit = MassInf[mask]
    errfit = MassErr[mask]
    models = []
    for signvec in itertools.product((-1,1), repeat=len(xfit)):
        x = np.array(list(signvec))
        Y = yfit + x*errfit
        # Y = np.zeros(len(xfit))
        # for i in range(len(x)):
            # Y[i] = yfit[i] -(-1)**(x[i])
        X = np.array([xfit]).transpose()
        models.append(LinearRegression().fit(X.reshape(-1,1), Y))

    ints = np.array([model.intercept_ for model in models])
    coefs = np.array([model.coef_[0] for model in models])

    def predict(x):
        return (x*np.mean(coefs) + np.mean(ints))

    gc = -np.sum(ints)/np.sum(coefs)

    gcs = -ints/coefs
    gcerr = (np.max(gc-gcs), np.max(gcs-gc))

    xlist = scipy.linspace(min(xfit), max(xfit), 100)
    axes[0].plot(xlist, predict(xlist), c='g')

    xlist = scipy.linspace(max(xfit), gc, 100)
    axes[0].plot(xlist, predict(xlist), c='g', linestyle='--')

    msg = [
            r"$g_c = {:.7f} +{:.7f} -{:.7f}$".format(gc, gcerr[1], gcerr[0])
    ]
    print(msg[0])

    for i, m in enumerate(msg):
        axes[0].text(0.8, 0.85-i*0.05, msg[i], horizontalalignment='center',
            verticalalignment='center', fontsize=13, transform=axes[0].transAxes)

    axes[1].errorbar(xfit, yfit-predict(xfit), errfit, color='r')
    axes[1].errorbar(xfit, MassInf2[mask]-predict(xfit), MassErr2[mask], color='b')


argv = sys.argv

if len(argv) < 1:
    print(argv[0], "")
    sys.exit(-1)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

f, axes = plt.subplots(2, 1, sharex='col')
f.subplots_adjust(hspace=0, wspace=0, top=0.94, right=0.95, left=0)

plotvsG(glist, axes)

print(gc)

axes[0].set_ylabel(r"$m_{ph}$")
axes[0].set_ylim(-0.01, 1)
axes[0].set_xlim(xmin, gc)
axes[1].set_ylabel("residuals")
axes[1].set_xlabel(r"$g$")

# Remove common tick
axes[1].yaxis.set_major_locator(MaxNLocator(prune='upper', nbins=8))

fname = ".{0}".format(output)

# MASS
plt.figure(1, figsize=(4, 2.5), dpi=300)
s = "MvsG"
if nparams==2:
    s += "_p=2"
plt.savefig(s+fname, bbox_inches='tight')
