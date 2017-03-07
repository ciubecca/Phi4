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

# x represent 1/ET^2
def fitfun1(x, a, b, c):
    return a + b*x**c

def fitfun2(x, a, b, c):
    return a + b*x + c*x**(3/2)

marker = 'o'
markersize = 2.5

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

# Weight function for linear regression
def sigma(x):
    return 1/x**2
sigmastr = "x-2"

# def sigma(x):
    # return 1
# sigmastr = "1"

# def sigma(x):
    # return 1/x
# sigmastr = "x-1"


def Massfun(L, a, b, c):
    return a + (b/L)*exp(-c*L)

def Lambdafun(L, a, b):
    return a - sqrt(b/(2*pi*L**3))*exp(-b*L)

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


Llist = [6, 7, 8, 9, 10]
ETs = {6:[10,27], 7:[10,25], 8:[10,22.5], 9:[10,20.5], 10:[10,19.5]}

ETlists = {L: scipy.linspace(ETs[L][0], ETs[L][1], (ETs[L][1]-ETs[L][0])*2+1)
        for L in Llist}

color = {6:"r", 7:"y", 8:"g", 9:"k", 10:"b"}

marker = 'o'
markersize = 2.5
ymin = [0, 0]
ymax = [0, 0]
xmax = max(Llist)+0.1
xmin = min(Llist)-0.1

db = database.Database()

def extrapolate(L):

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
    data = spectrum[1][ren][:,0]/L
    popt, pcov = curve_fit(fitfun, 1/xlist**2, data, sigma=sigma(xlist),
            absolute_sigma=absolute_sigma)
    err = errcoeff*np.sqrt(pcov[0,0])
    LambdaBar = (popt[0], err)


    # MASS
    data = masses[-1][ren][0]
    popt, pcov = curve_fit(fitfun, 1/xlist**2, data, sigma=sigma(xlist), p0=[0,0,1],
            absolute_sigma=absolute_sigma)
    err = errcoeff*np.sqrt(pcov[0,0])
    MassBar = (popt[0], err)

    return LambdaBar, MassBar


def plotvsL(Llist):
    xlist = Llist

    LambdaVec = scipy.zeros(len(Llist))
    LambdaErrVec = scipy.zeros(len(Llist))
    MassVec = scipy.zeros(len(Llist))
    MassErrVec = scipy.zeros(len(Llist))

    for i, L in enumerate(Llist):
        (Lambda, LambdaErr), (Mass, MassErr) = extrapolate(L)
        LambdaVec[i] = Lambda
        MassVec[i] = Mass
        LambdaErrVec[i] = LambdaErr
        MassErrVec[i] = MassErr

    ymax[0] = max(LambdaVec+LambdaErrVec)
    ymin[0] = min(LambdaVec-LambdaErrVec)
    ymax[1] = max(MassVec+MassErrVec)
    ymin[1] = min(MassVec-MassErrVec)


    plt.figure(1)
    plt.errorbar(xlist, LambdaVec, yerr=LambdaErrVec, linestyle=' ', marker=marker)
    # Fit with exponential
    popt, pcov = curve_fit(Lambdafun, xlist, LambdaVec, sigma=LambdaErrVec,
            absolute_sigma=absolute_sigma, p0=[-1, MassVec[-1]])
    print(popt)
    print("Best fit for m_ph from Lambda fit: {} +- {}", popt[1], np.sqrt(pcov[1,1]))
    xdata = scipy.linspace(xmin, xmax, 100)
    plt.plot(xdata, Lambdafun(xdata, *popt))

    plt.figure(2)
    plt.errorbar(xlist, MassVec, yerr=MassErrVec, linestyle=' ', marker=marker)
    # Fit with exponential
    return
    popt, pcov = curve_fit(Massfun, xlist, MassVec, sigma=MassErrVec,
            absolute_sigma=absolute_sigma, p0=[MassVec[-1],-1,1.5])
    print(popt)
    print("Best fit for m_ph from Mass fit: {} +- {}", popt[0], np.sqrt(pcov[0,0]))
    xdata = scipy.linspace(xmin, xmax, 100)
    plt.plot(xdata, Massfun(xdata, *popt))

argv = sys.argv

if len(argv) < 2:
    print(argv[0], "<g>")
    sys.exit(-1)

g = float(argv[1])

print("g=", g)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))

plotvsL(Llist)

title = r"$g$={0:.1f}, {1}, $\sigma$={2}".format(g, fstr(), sigmastr)
fname = "g={0:.1f}_ftype={1}_sigma={2}.{3}".format(g,ftype,sigmastr,output)
loc = "upper left"


# VACUUM ENERGY DENSITY
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$E_0/L$")
ymargin = 10**(-5)
plt.xlim(xmin, xmax)
plt.ylim(ymin[0]-ymargin, ymax[0]+ymargin)
# plt.legend(loc=loc)
plt.savefig("extrLambdavsL_"+fname)


# MASS
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$M$")
plt.xlim(xmin, xmax)
plt.ylim(ymin[1]-ymargin, ymax[1]+ymargin)
# plt.legend(loc=loc)
plt.savefig("extrMassvsL_"+fname)
