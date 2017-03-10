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
from extrapolate import Extrapolator, ETmax

Llist = [5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]


mph = {k:None for k in (-1,1)}
mpherr = {k: None for k in (-1,1)}

marker = 'o'
markersize = 2.5

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))

def Massfun(L, a, b, c):
    return a + (b/L)*exp(-c*a*sqrt(3/2)*L)
boundsMass = ([0,0,1],[1,np.inf,np.inf])

def Lambdafun(L, a, b):
    return a - sqrt(sqrt(b**2)/(2*pi*L**3))*exp(-b*L)
boundsLambda = ([-np.inf,0],[0,1])

output = "png"
renlist = ("rentails", "renloc")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

klist = (1,-1)

ymin = {}
ymax = {}
xmax = max(Llist)+0.1
xmin = min(Llist)-0.1

db = database.Database()


def plotvsL(Llist):
    xlist = np.array(Llist)
    global mph
    global mpherr

    Lambda = {ren: np.zeros(len(Llist)) for ren in renlist}
    LambdaInf = np.zeros(len(Llist))
    Mass = {ren: np.zeros(len(Llist)) for ren in renlist}
    MassInf = np.zeros(len(Llist))

    for i,L in enumerate(Llist):
        for ren in renlist:
            E0 = db.getEigs(1, ren, g, L, ETmax[L])[0]
            E1 = db.getEigs(-1, ren, g, L, ETmax[L])[0]
            Lambda[ren][i] = E0/L
            Mass[ren][i] = (E1-E0)

        e = {}
        e[1] = Extrapolator(db, 1, L, g)
        e[-1] = Extrapolator(db, -1, L, g)
        e[1].train(alpha)
        e[-1].train(alpha)
        LambdaInf[i] = e[1].asymptoticValue()/L
        MassInf[i] = e[-1].asymptoticValue()-e[1].asymptoticValue()

    ymax[1] = max(max(max(Lambda[ren]) for ren in renlist), max(LambdaInf))
    ymin[1] = min(min(min(Lambda[ren]) for ren in renlist), min(LambdaInf))
    ymax[-1] = max(max(max(Mass[ren]) for ren in renlist), max(MassInf))
    ymin[-1] = min(min(min(Mass[ren]) for ren in renlist), min(MassInf))

    # print(",".join(map(str,list(Llist))))
    # print(",".join(map(str,list(LambdaInf))))
    # print(",".join(map(str,list(MassInf))))

    # Lambda
    plt.figure(1)
    # Plot non-extrapolated data
    for ren in renlist:
        plt.plot(xlist, Lambda[ren], marker=marker, label=ren)

    # Plot extrapolated data
    plt.scatter(xlist, LambdaInf, marker=marker, label=r"$E_T=\infty$")

    err = np.array([0.01]*len(Llist))
    popt, pcov = curve_fit(Lambdafun, xlist.ravel(), LambdaInf.ravel(),
            bounds=boundsLambda, method='dogbox')
    print("Best fit values for Lambda fit:", popt)
    mph[1] = popt[1]
    mpherr[1] = np.sqrt(pcov[1,1])
    xdata = scipy.linspace(xmin, xmax, 100)
    plt.plot(xdata, Lambdafun(xdata, *popt))

    # Mass
    plt.figure(2)
    # Plot non-extrapolated data
    for ren in renlist:
        plt.plot(xlist, Mass[ren], marker=marker, label=ren)

    # Plot extrapolated data
    plt.scatter(xlist, MassInf, marker=marker, label=r"$E_T=\infty$")

    popt, pcov = curve_fit(Massfun, xlist.ravel(), MassInf.ravel(),
            bounds=boundsMass, method='dogbox')
    print("Best fit values for Mass fit:", popt)
    mph[-1] = popt[0]
    mpherr[-1] = np.sqrt(pcov[0,0])
    xdata = scipy.linspace(xmin, xmax, 100)
    plt.plot(xdata, Massfun(xdata, *popt))


argv = sys.argv

if len(argv) < 3:
    print(argv[0], "<g> <alpha>")
    sys.exit(-1)

g = float(argv[1])
alpha = float(argv[2])

print("g=", g)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plotvsL(Llist)

fname = "g={0:.1f}_alpha={1}.{2}".format(g, alpha, output)
loc = "upper left"


# VACUUM ENERGY DENSITY
title = r"$g$={:.1f}, $\alpha$={}, $m_{{ph}}={}\pm{}$".format(g, alpha,
        mph[1], mpherr[1])
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$E_0/L$")
ymargin = 10**(-5)
plt.xlim(xmin, xmax)
plt.ylim(ymin[1]-ymargin, ymax[1]+ymargin)
plt.legend(loc=loc)
plt.savefig("extrLambdavsL_"+fname)

# MASS
title = r"$g$={:.1f}, $\alpha$={}, $m_{{ph}}={}\pm{}$".format(g, alpha,
        mph[-1], mpherr[-1])
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$M$")
plt.xlim(xmin, xmax)
plt.ylim(ymin[-1]-ymargin, ymax[-1]+ymargin)
plt.legend(loc=loc)
plt.savefig("extrMassvsL_"+fname)
