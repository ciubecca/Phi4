import sys
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
from cycler import cycler
import database
from sys import exit
import numpy as np
from extrapolate import Extrapolator

xmargin = 10**(-4)

power = {"renloc":2, "rentails":3}

Llist = [6, 8, 10]

output = "pdf"
renlist = ("rentails",)

plt.style.use('ggplot')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
params = {'legend.fontsize': 8}
plt.rcParams.update(params)


def featureVec(ET):
    return [1/ET**3]

neigs = 6

color = {5:"k", 6:"r", 6.5:"k", 7:"g", 8:"y", 9:"b", 10:"c"}

marker = 'o'
ymax = {1:-10, -1:0}
ymin = {1:10, -1:10}
xmax = {"rentails":0, "renloc":0}

db = database.Database("data/spectra3.db")


def plotvsET(Llist, axes, g, alpha):

    global ymin, ymax, xmax

    for ren in renlist:
        for i,L in enumerate(Llist):
            spectrum = {}
            res = {}

            for k in (-1,1):
                e = Extrapolator(db, k, L, g, ren=ren)
                ETlist = e.ETlist
                xlist = 1/ETlist**power[ren]
                xmax[ren] = max(max(xlist), xmax[ren])
                spectrum[k] = e.spectrum

                e.train(alpha=alpha, featureVec=featureVec)
                res[k] = spectrum[k] - e.predict(xlist**(-1/power[ren]))


            label = "ren = {}, L = {}".format(ren,L)

            " VACUUM ENERGY RESIDUALS"
            ax = axes[0]
            ax.plot(1/xlist**(1/3),
                    res[1]/L, color=color[L], label=label, marker='o',
                    markersize=5)

            ymax[1] = max(ymax[1], max(res[1]/L))
            ymin[1] = min(ymin[1], min(res[1]/L))

            " MASS RESIDUALS"
            ax = axes[1]
            ax.plot(1/xlist**(1/3), res[1]-res[-1], color=color[L], label=label,
                    marker='o', markersize=5)

            ymax[-1] = max(ymax[-1], max(res[1]-res[-1]))
            ymin[-1] = min(ymin[-1], min(res[1]-res[-1]))



def main(argv):

    argv = sys.argv

    if len(argv) < 3:
        print(argv[0], "<g> <alpha>")
        sys.exit(-1)

    g = float(argv[1])
    alpha = float(argv[2])

    print("g=", g)


    title = r"$g$={0:.1f}, $\alpha$={1}".format(g, alpha)
    fname = "g={0:.1f}_alpha={1}.{2}".format(g,alpha,output)

    f, axes = plt.subplots(2, 1, sharex='col')
    f.subplots_adjust(hspace=0, wspace=0, top=0.93, right=0.95, left=0.1)
    f.suptitle(r"$g={}, \quad \alpha={}$".format(g,alpha), fontsize=15)

    plotvsET(Llist, axes, g, alpha)


    axes[0].legend(loc=2)
    # axes[0].set_xlim(0, xmax["rentails"]+10**(-5))
    axes[0].legend(loc=1)
    axes[0].set_ylabel(r"$E_0/L$")
    ymargin = (ymax[1]-ymin[1])/100
    axes[0].set_ylim(ymin[1]-ymargin, ymax[1]+ymargin)
    axes[0].invert_xaxis()

    axes[1].set_ylim(ymin[-1]-ymargin, ymax[-1]+ymargin)
    axes[1].set_xlabel(r"$1/E_{{T}}^{}$".format(power["rentails"]))
    axes[1].set_ylabel(r"$E_1-E_0$")


# JOINT PLOT
    plt.savefig("resvsET_"+fname)


if __name__=="__main__":
    main(sys.argv)
