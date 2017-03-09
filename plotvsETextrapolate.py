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

alpha = {1:0.1, 2:0.01}

Llist = [6, 7, 8, 9, 10]

output = "png"
renlist = ("raw", "rentails", "renloc")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

neigs = 6

color = {6:"r", 7:"y", 8:"g", 9:"k", 10:"b"}

marker = 'o'
markersize = 2.5
ymax = {1:-10, -1:0}
ymin = {1:10, -1:10}
xmax = 0

db = database.Database()

fignum= {1:1, -1:2}

def plotvsET(L):

    spectrum = {}
    ydata = {}
    yinf = {}

    global ymin, ymax, xmax

    for k in (-1,1):
        e = Extrapolator(db, k, "rentails", L, g)
        ETlist = e.ETlist
        xlist = 1/ETlist**2
        xmax = max(max(xlist), xmax)
        spectrum[k] = e.spectrum

        e.train(alpha=alpha[g])
        xdata = scipy.linspace(0, min(ETlist)**-2, 100)
        ydata[k] = e.predict(xdata**(-1/2))

        yinf[k] = e.asymptoticValue()


    label = "L = {}".format(L)
    " VACUUM ENERGY "
    plt.figure(fignum[1])
    plt.scatter(xlist, spectrum[1]/L, label=label, marker=marker, color=color[L])
    plt.plot(xdata, ydata[1]/L, color=color[L])
    ymax[1] = max(ymax[1], max(ydata[1])/L)
    ymin[1] = min(ymin[1], min(ydata[1])/L)

    " MASS "
    plt.figure(fignum[-1])
    plt.scatter(xlist, spectrum[-1]-spectrum[1], label=label,
            marker=marker, color=color[L])
    plt.plot(xdata, ydata[-1]-ydata[1], color=color[L])
    ymax[-1] = max(ymax[-1], max(ydata[-1]-ydata[1]))
    ymin[-1] = min(ymin[-1], min(ydata[-1]-ydata[1]))


argv = sys.argv

if len(argv) < 2:
    print(argv[0], "<g>")
    sys.exit(-1)

g = float(argv[1])

print("g=", g)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)


for L in Llist:
    plotvsET(L)


title = r"$g$={0:.1f}, $\alpha$={1}".format(g, alpha[g])
fname = "g={0:.1f}_alpha={1}.{2}".format(g,alpha[g],output)
loc = "upper left"


# VACUUM ENERGY DENSITY
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$1/E_{T}^2$")
plt.ylabel(r"$E_0/L$")
plt.xlim(0-xmargin, xmax+xmargin)
ymargin = (ymax[1]-ymin[1])/100
plt.ylim(ymin[1]-ymargin, ymax[1]+ymargin)
plt.legend(loc=loc)
plt.savefig("extrLambda_"+fname)

# M
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$1/E_{T}^2$")
plt.ylabel(r"$M$")
plt.xlim(0-xmargin, xmax+xmargin)
ymargin = (ymax[-1]-ymin[-1])/100
plt.ylim(ymin[-1]-ymargin, ymax[-1]+ymargin)
plt.legend(loc=loc)
plt.savefig("extrMass"+fname)
