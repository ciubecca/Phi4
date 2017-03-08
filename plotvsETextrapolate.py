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
from sklearn.linear_model import Ridge


# x represent 1/ET^2
def fitfun1(x, a, b, c):
    return a + b*x**c

def fitfun2(x, a, b, c):
    return a + b*x + c*x**(3/2)
def fitfun3(x, a, b):
    return a + b*x

xmargin = 10**(-4)
ymargin = 10**(-6)

alpha = {1:0.1, 2:0.01}

Llist = [6, 7, 8, 9, 10]
ETs = {5:30.5, 5.5:28.5, 6:28, 6.5:24.5, 7:25, 7.5:22, 8:23,
        8.5:20, 9:21, 9.5:20, 10:20}

Erange = 10

absolute_sigma=False
errcoeff = 3.

ftype = 2
if ftype==1:
    fitfun = fitfun1
    def fstr():
        return "f(x) = {} + {}*x**{}".format("a","b","c")
elif ftype==2:
    fitfun = fitfun2
    def fstr():
        return "f(x) = {} + {}*x +{}*x**(3/2)".format("a","b","c")
elif ftype==3:
    fitfun = fitfun3
    def fstr():
        return "f(x) = {} + {}*x".format("a","b")

# Weight function for linear regression
def weights(ET):
    # return abs(1/ET**2)
    return 1
sigmastr = "1"


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



ETlists = {L: scipy.linspace(ETs[L]-Erange, ETs[L], (Erange)*2+1)
        for L in Llist}

color = {6:"r", 7:"y", 8:"g", 9:"k", 10:"b"}

marker = 'o'
markersize = 2.5
ymin = [0, 10]
ymax = [-10, 0]
xmax = max(max(1/ETlists[L]**2) for L in Llist)

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

    popt, pcov = curve_fit(fitfun, 1/xlist**2, data, sigma=weights(xlist),
            absolute_sigma=absolute_sigma)
    print("L={}, Lambda, popt={}".format(L,popt))
    err = errcoeff*np.sqrt(pcov[0,0])

    X = scipy.array([1/xlist**2, 1/xlist**(3)]).transpose()
    model = Ridge(alpha=alpha[g], normalize=True)
    model.fit(X, data, sample_weight=1/weights(xlist))

    xdata = scipy.linspace(0, xmax, 100)
    ydata = model.predict(scipy.array([xdata, xdata**(3/2)]).transpose())
    # ydata = fitfun(xdata, *popt)
    plt.plot(xdata, ydata, color=color[L])
    # plt.errorbar([0], [popt[0]], yerr=err, color=color[L])
    ymax[0] = max(ymax[0], max(data))
    # ymin[0] = min(ymin[0], min(data), popt[0]-err)
    ymin[0] = min(ymin[0], min(data), popt[0])

    # MASS
    plt.figure(2)
    data = masses[-1][ren][0]
    label = "L = {}".format(L)
    plt.scatter(1/xlist**2, data, label=label, marker=marker, color=color[L])

    # popt, pcov = curve_fit(fitfun, 1/xlist**2, data, sigma=sigma(xlist), p0=[0,0,1],
            # absolute_sigma=absolute_sigma)
    # print("L={}, Mass, popt={}".format(L,popt))
    # err = errcoeff*np.sqrt(pcov[0,0])
    # xdata = scipy.linspace(0, xmax, 100)
    # plt.plot(xdata, fitfun(xdata, *popt), color=color[L])
    # plt.errorbar([0], [popt[0]], yerr=err, color=color[L])
    # ymax[1] = max(ymax[1], max(data))
    # ymin[1] = min(ymin[1], min(data), popt[0]-err)


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


title = r"$g$={0:.1f}, $\alpha$={1}".format(g, alpha[g])
fname = "g={0:.1f}_alpha={1}.{2}".format(g,alpha[g],output)
loc = "upper left"


# VACUUM ENERGY DENSITY
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$1/E_{T}^2$")
plt.ylabel(r"$E_0/L$")
plt.xlim(0-xmargin, xmax+xmargin)
plt.ylim(ymin[0]-ymargin, ymax[0]+ymargin)
plt.legend(loc=loc)
plt.savefig("extrLambda_"+fname)

sys.exit(0)

# MASS
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$1/E_{T}^2$")
plt.ylabel(r"$M$")
plt.xlim(0-xmargin, xmax+xmargin)
plt.ylim(ymin[1]-ymargin, ymax[1]+ymargin)
plt.legend(loc=loc)
plt.savefig("extrMass_"+fname)
