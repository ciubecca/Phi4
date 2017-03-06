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

# x represent 1/ET^2
def fitfun1(x, a, b, c):
    return a + b*x**c

def fitfun2(x, a, b, c):
    return a + b*x + c*x**(3/2)

absolute_sigma=False

ftype = 2
if ftype==1:
    fitfun = fitfun1
    def fstr():
        return "f(x) = {} + {}*x**{}".format("a","b","c")
elif ftype==2:
    fitfun = fitfun2
    def fstr():
        return "f(x) = {} + {}*x +{}*x**(3/2)".format("a","b","c")

# Weight function for linear regression
def sigma(x):
    return 1/x**2
sigmastr = "x-2"

def sigma(x):
    return 1
sigmastr = "1"

def sigma(x):
    return 1/x
sigmastr = "x-1"

output = "png"
renlist = ("raw", "rentails", "renloc")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

klist = (1,-1)

neigs = 6

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 1.5


Llist = [7, 9, 10]
ETs = {6:[10,26.5], 7:[10,24.5], 8:[10,22.5], 9:[10,20.5], 10:[10,19.5]}

ETlists = {L: scipy.linspace(ETs[L][0], ETs[L][1], (ETs[L][1]-ETs[L][0])*2+1)
        for L in Llist}

color = {6:"r", 7:"y", 8:"g", 9:"k", 10:"b"}

marker = 'o'
markersize = 2.5
ymin = [0, 10]
ymax = [-10, 0]
xmax = max(1/ETs[L][0]**2 for L in Llist)

db = database.Database()

def plotvsET(L):

    global ymin, ymax

    xlist = ETlists[L]

    spectrum = {k:{ren:[] for ren in renlist} for k in klist}
    masses = {k:{ren:[] for ren in renlist} for k in klist}

    ren = "rentails"

    for k in klist:
        for ET in xlist:
            EL = ratioELET*ET
            ELp = ratioELpET*ET
            ELpp = ratioELppELp*ELp

            approxQuery = {"g":g, "L":L, "ET":ET}
            exactQuery = {"k": k, "ren":ren, "neigs":neigs}
            boundQuery = {}

            if ren=="rentails":
                exactQuery["maxntails"] = None
                exactQuery["tailsComputedAtET"] = ET
                approxQuery["EL"] = EL
                approxQuery["ELp"] = ELp
                approxQuery["ELpp"] = ELpp

            try:
                spectrum[k][ren].append(db.getObjList('spec', exactQuery,
                    approxQuery, boundQuery, orderBy="date")[0])

            except IndexError:
                print("Not found:", exactQuery, approxQuery)
                exit(-1)

    # Convert to array
    for k in klist:
        spectrum[k][ren] = array(spectrum[k][ren])

    # Mass
    for k in (-1,1):
        for i in range(int((1+k)/2), neigs):
            masses[k][ren].append(spectrum[k][ren][:,i]-spectrum[1][ren][:,0])

    # Convert to array
    for k in klist:
        masses[k][ren] = array(masses[k][ren])


    # VACUUM ENERGY
    plt.figure(1)
    data = spectrum[1][ren][:,0]/L
    label = "L = {}".format(L)
    plt.scatter(1/xlist**2, data, label=label, marker=marker,
            color=color[L])

    popt, pcov = curve_fit(fitfun, 1/xlist**2, data, sigma=sigma(xlist),
            absolute_sigma=absolute_sigma)
    # err = max(abs(fitfun(x,*popt)-y) for x,y in list(zip(1/xlist**2.,data))[-10:])
    err = 2*np.sqrt(pcov[0,0])
    xdata = scipy.linspace(0, xmax, 100)
    plt.plot(xdata, fitfun(xdata, *popt), color=color[L])
    plt.errorbar([0], [popt[0]], yerr=err, color=color[L])
    ymax[0] = max(ymax[0], max(data))
    ymin[0] = min(ymin[0], min(data), popt[0]-err)

    # MASS
    plt.figure(2)
    data = masses[-1][ren][0]
    label = "L = {}".format(L)
    plt.scatter(1/xlist**2, data, label=label, marker=marker, color=color[L])

    popt, pcov = curve_fit(fitfun, 1/xlist**2, data, sigma=sigma(xlist), p0=[0,0,1],
            absolute_sigma=absolute_sigma)
    print(popt)
    # err = max(abs(fitfun(x,*popt)-y) for x,y in list(zip(1/xlist**2.,data))[-10:])
    err = 2*np.sqrt(pcov[0,0])
    xdata = scipy.linspace(0, xmax, 100)
    plt.plot(xdata, fitfun(xdata, *popt), color=color[L])
    plt.errorbar([0], [popt[0]], yerr=err, color=color[L])
    ymax[1] = max(ymax[1], max(data))
    ymin[1] = min(ymin[1], min(data), popt[0]-err)


argv = sys.argv

if len(argv) < 2:
    print(argv[0], "<g>")
    sys.exit(-1)

g = float(argv[1])

print("g=", g)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))

for L in Llist:
    plotvsET(L)


title = r"$g$={0:.1f}, {1}, $\sigma$={2}".format(g, fstr(), sigmastr)
fname = "g={0:.1f}_ftype={1}_sigma={2}.{3}".format(g,ftype,sigmastr,output)
loc = "upper left"


# VACUUM ENERGY DENSITY
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$1/E_{T}^2$")
plt.ylabel(r"$E_0/L$")
margin = 10**(-4)
plt.xlim(0-margin, xmax+margin)
plt.ylim(ymin[0]-margin, ymax[0]+margin)
plt.legend(loc=loc)
plt.savefig("extrLambda_"+fname)


# MASS
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$1/E_{T}^2$")
plt.ylabel(r"$M$")
plt.xlim(0-margin, xmax+margin)
plt.ylim(ymin[1]-margin, ymax[1]+margin)
plt.legend(loc=loc)
plt.savefig("extrMass_"+fname)
