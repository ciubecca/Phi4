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
from matplotlib.ticker import MaxNLocator
from itertools import cycle

power = {"raw":2, "renloc":2, "rentails":3}
labelstr = {"raw":"raw", "renloc":"local LO", "rentails":"NLO"}

Llist = [6, 8, 10]

output = "pdf"
renlist = ("raw", "rentails", "renloc")

plt.style.use('ggplot')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
params = {'legend.fontsize': 8}
plt.rcParams.update(params)

# plt.rc('axes', prop_cycle=(
    # cycler('marker', ['x', 'o', 'v'])
    # +cycler('linestyle', ['-', '--', ':'])
    # +cycler('markersize', [5.,3.,3.])
    # +cycler('color', ['r','b','g'])
    # ))

sty_cycle = cycler('marker', ['D', 'o', '^'])\
    +cycler('linestyle', ['-', '--', ':'])\
    +cycler('s', [10.,8.,20.])\
    +cycler('color', ['r','b','g'])

sty_cycle_loc = cycler('marker', ['x', '*', 's'])\
    +cycler('linestyle', ['-', '--', ':'])\
    +cycler('s', [12.,20.,15.])\
    +cycler('color', ['k','m','c'])


neigs = 6

color = {5:"k", 6:"r", 6.5:"k", 7:"g", 8:"y", 9:"b", 10:"c"}

marker = 'o'
ymax = {1:-10, -1:0}
ymin = {1:10, -1:10}
xmax = {"rentails":0, "renloc":0, "raw":0}

db = database.Database("data/spectra3.db")


def plotvsET(Llist, axes):

    global ymin, ymax, xmax

    for ren in renlist:

        if ren=="rentails":
            j=0
        else:
            j=1

        if ren=="renloc":
            optIter = cycle(sty_cycle_loc)
        else:
            optIter = cycle(sty_cycle)

        for i,L in enumerate(Llist):

            opt = next(optIter)

            spectrum = {}
            ydata = {}
            yerr = {}
            yinf = {}

            e = Extrapolator(db, L, g, ren=ren)
            e.train(neigs=1)
            ETlist = e.ETlist
            xlist = 1/ETlist**power[ren]
            xmax[ren] = max(max(xlist), xmax[ren])

            for k in (-1,1):
                spectrum[k] = e.spectrumSub[k][:,0]

                if ren=="rentails":
                    xdata = scipy.linspace(0, min(ETlist)**-power[ren], 100)
                    ydata[k] = e.predict(k, xdata**(-1/power[ren]))[0]
                    yinf[k] = e.asymValue(k)[0]
                    yerr[k] = e.asymErr(k)[0]

            # Position of error bar
            xerr = np.array([10**(-5)])
            xerr = np.array([0])

            label = "{}, L = {}".format(labelstr[ren],L)
            " VACUUM ENERGY "
            ax = axes[j]
            ax.scatter(xlist, spectrum[1]/L, label=label, **opt)
            ymax[1] = max(ymax[1], max(spectrum[1]/L))
            ymin[1] = min(ymin[1], min(spectrum[1]/L))
            if ren=="rentails":
                ax.plot(xdata, ydata[1]/L, c=opt['color'])

                err = np.expand_dims(yerr[1], axis=1)
                # ax.errorbar(xerr, np.array([yinf[1]/L]), err/L, color=color[L],
                    # elinewidth=1)

                ymax[1] = max(ymax[1], max(ydata[1]+err[0])/L, max(ydata[1])/L)
                ymin[1] = min(ymin[1], min(ydata[1]-err[1])/L, min(ydata[1])/L)


argv = sys.argv

if len(argv) < 2:
    print(argv[0], "<g>")
    sys.exit(-1)

g = float(argv[1])

print("g=", g)


title = r"$g$={0:.1f}".format(g)
fname = "g={0:.1f}.{1}".format(g,output)


f, axes = plt.subplots(1, 2, sharey='row')
f.subplots_adjust(hspace=0, wspace=0, top=.94, right=1, left=0)
f.suptitle(r"$g={}$".format(g), fontsize=15)

plotvsET(Llist, axes)


axes[0].set_xlim(0, xmax["rentails"]+10**(-4))
axes[0].legend(loc=1, fontsize=10)
axes[1].set_xlim(0, xmax["renloc"]+10**(-5))
axes[1].legend(loc=2, fontsize=10)
axes[0].set_ylabel(r"$\mathcal{E}_0/L$", fontsize=12)
ymargin = (ymax[1]-ymin[1])/100
axes[0].set_ylim(ymin[1]-ymargin, ymax[1]+ymargin)
axes[0].invert_xaxis()

axes[0].set_xlabel(r"$1/E_{{T}}^{}$".format(power["rentails"]), fontsize=15)
axes[1].set_xlabel(r"$1/E_{{T}}^{}$".format(power["renloc"]), fontsize=15)

# Remove common tick
axes[1].xaxis.set_major_locator(MaxNLocator(prune='lower', nbins=6))

# JOINT PLOT
plt.savefig("fitvsETRawE0"+fname, bbox_inches='tight')
