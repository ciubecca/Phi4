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

power = 3

Llist = [6, 6.5, 7, 8, 9, 10]

output = "pdf"
renlist = ("raw", "rentails", "renloc")

plt.style.use('ggplot')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
params = {'legend.fontsize': 8}
plt.rcParams.update(params)

neigs = 6

color = {5:"k", 6:"r", 6.5:"k", 7:"g", 8:"y", 9:"b", 10:"c"}

marker = 'o'
ymax = {1:-10, -1:0}
ymin = {1:10, -1:10}
xmax = 0

db = database.Database()


def plotvsET(Llist, axes):

    for i,L in enumerate(Llist):
        spectrum = {}
        ydata = {}
        yinf = {}

        global ymin, ymax, xmax

        for k in (-1,1):
            e = Extrapolator(db, k, L, g)
            ETlist = e.ETlist
            xlist = 1/ETlist**power
            xmax = max(max(xlist), xmax)
            spectrum[k] = e.spectrum

            e.train(alpha=alpha)
            xdata = scipy.linspace(0, min(ETlist)**-power, 100)
            ydata[k] = e.predict(xdata**(-1/power))

            yinf[k] = e.asymValue()


        label = "L = {}".format(L)
        " VACUUM ENERGY "
        if i%2==0:
            ax = axes[0,0]
        else:
            ax = axes[0,1]
        ax.scatter(xlist, spectrum[1]/L, label=label, marker=marker, color=color[L])
        ax.plot(xdata, ydata[1]/L, color=color[L])
        ymax[1] = max(ymax[1], max(ydata[1])/L, max(spectrum[1]/L))
        ymin[1] = min(ymin[1], min(ydata[1])/L, min(spectrum[1]/L))

        " MASS "
        if i%2==0:
            ax = axes[1,0]
        else:
            ax = axes[1,1]
        mass = spectrum[-1]-spectrum[1]
        ax.scatter(xlist, mass, label=label, marker=marker, color=color[L])
        ax.plot(xdata, ydata[-1]-ydata[1], color=color[L])
        ymax[-1] = max(ymax[-1], max(ydata[-1]-ydata[1]), max(mass))
        ymin[-1] = min(ymin[-1], min(ydata[-1]-ydata[1]), min(mass))


argv = sys.argv

if len(argv) < 3:
    print(argv[0], "<g> <alpha>")
    sys.exit(-1)

g = float(argv[1])
alpha = float(argv[2])

print("g=", g)


title = r"$g$={0:.1f}, $\alpha$={1}".format(g, alpha)
fname = "g={0:.1f}.{1}".format(g,output)


f, axes = plt.subplots(2, 2, sharex='col', sharey='row')
f.subplots_adjust(hspace=0, wspace=0, top=0.93, right=0.95)
f.suptitle(r"$g={}, \quad \alpha={}$".format(g,alpha), fontsize=15)

plotvsET(Llist, axes)


axes[0,1].set_xlim(0, xmax+xmargin)
axes[0,1].legend(loc=2)
axes[0,0].set_xlim(0, xmax+xmargin)
axes[0,0].legend(loc=1)
axes[0,0].set_ylabel(r"$E_0/L$")
ymargin = (ymax[1]-ymin[1])/100
axes[0,0].set_ylim(ymin[1]-ymargin, ymax[1]+ymargin)
axes[0,0].invert_xaxis()

axes[1,0].set_ylim(ymin[-1]-ymargin, ymax[-1]+ymargin)
axes[1,0].set_xlabel(r"$1/E_{{T}}^{}$".format(power))
axes[1,1].set_xlabel(r"$1/E_{{T}}^{}$".format(power))
axes[1,0].set_ylabel(r"$E_1-E_0$")


# JOINT PLOT
plt.savefig("fitvsET_"+fname)
