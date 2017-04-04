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
from scipy.special import kn

Llist = {}
Llist["rentails"] = [5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]
Llist["renloc"] = [6,8,10]

output = "pdf"
renlist = ("rentails", "renloc")
renlist = ("rentails", )

method='dogbox'

marker = 'o'
markersize = 2.5

plt.style.use('ggplot')
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))

# def Massfun(L, a, b):
    # return a - (3*b/(16*pi**2*L))*exp(-a*L)
# boundsMass = ([0,-np.inf],[1,np.inf])
# fmassStr = r"$m_{ph} - B \frac{3}{16 \pi^2 L} e^{- m_{ph} L}$"


def Massfun(L, m, b, c):
    return m - m*1/(L*m)**(3/2)*exp(-m*L)*(b+c/(L*m))
boundsMass = ([0,-np.inf,-np.inf],[1,np.inf,np.inf])
fmassStr = r"$m_{ph} - \frac{m_{ph}}{(L m_{ph})^{3/2}} e^{- m_{ph} L}(b+c/Lm)$"

def Lambdafun(L, a, b):
    return a - b/(pi*L)*kn(1, b*L)
boundsLambda = ([-np.inf,0],[0,1])
p0Lambda = [-.182, 0.5]
p0Lambda = [-.007, 0.93]
fvacStr = r"$\Lambda - \frac{m_{ph}}{\pi L} K_1(m_{ph} L)$"


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

klist = (1,-1)

ymin = {}
ymax = {}
xmax = max(Llist["rentails"])+0.1
xmin = min(Llist["rentails"])-0.1

db = database.Database("data/spectra3.db")


def plotvsL(Llist):

    Lambda = {ren: np.zeros(len(Llist[ren])) for ren in renlist}
    LambdaInf = np.zeros(len(Llist["rentails"]))
    LambdaErr = np.zeros(len(Llist["rentails"]))
    Mass = {ren: np.zeros(len(Llist[ren])) for ren in renlist}
    MassInf = np.zeros(len(Llist["rentails"]))
    MassErr = np.zeros(len(Llist["rentails"]))

    for ren in renlist:
        for i,L in enumerate(Llist[ren]):
            E0 = db.getEigs(1, ren, g, L, ETmax[ren][L])[0]
            E1 = db.getEigs(-1, ren, g, L, ETmax[ren][L])[0]
            Lambda[ren][i] = E0/L
            Mass[ren][i] = (E1-E0)

    for i,L in enumerate(Llist["rentails"]):
        e = {}
        e[1] = Extrapolator(db, 1, L, g)
        e[-1] = Extrapolator(db, -1, L, g)
        e[1].train(alpha)
        e[-1].train(alpha)
        LambdaInf[i] = e[1].asymValue()/L
        LambdaErr[i] = e[1].asymErr()/L
        MassInf[i] = e[-1].asymValue()-e[1].asymValue()
        MassErr[i] = max(e[-1].asymErr(),e[1].asymErr())

    ymax[1] = max(max(max(Lambda[ren]) for ren in renlist), max(LambdaInf+LambdaErr))
    ymin[1] = min(min(min(Lambda[ren]) for ren in renlist), min(LambdaInf-LambdaErr))
    ymax[-1] = max(max(max(Mass[ren]) for ren in renlist), max(MassInf+MassErr))
    ymin[-1] = min(min(min(Mass[ren]) for ren in renlist), min(MassInf-MassErr))


    # Lambda
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    # Plot non-extrapolated data
    for ren in renlist:
        plt.plot(Llist[ren], Lambda[ren], marker=marker, label=ren)

    # Plot extrapolated data
    # ax.scatter(xlist, LambdaInf, marker=marker, label=r"$E_T=\infty$")
    ax.errorbar(Llist["rentails"],
            LambdaInf, LambdaErr, marker=marker, label=r"$E_T=\infty$")

    # sigma = LambdaErr
    sigma = None
    popt, pcov = curve_fit(Lambdafun, Llist["rentails"], LambdaInf.ravel(),
            bounds=boundsLambda, method=method, p0=p0Lambda, sigma=sigma)
    xdata = scipy.linspace(xmin, xmax, 100)
    ax.plot(xdata, Lambdafun(xdata, *popt))
    msg = [
            r"$\Lambda = {:.7f} \pm {:.7f}$".format(popt[0],np.sqrt(pcov[0,0])),
            r"$m_{{ph}} = {:.7f} \pm {:.7f}$".format(popt[1],np.sqrt(pcov[1,1]))
            ]
    for i, m in enumerate(msg):
        ax.text(0.8, 0.2-i*0.05, msg[i], horizontalalignment='center',
            verticalalignment='center', fontsize=13, transform=ax.transAxes)


    # Mass
    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    # Plot non-extrapolated data
    for ren in renlist:
        plt.plot(Llist[ren], Mass[ren], marker=marker, label=ren)

    # Plot extrapolated data
    # ax.scatter(xlist, MassInf, marker=marker, label=r"$E_T=\infty$")
    ax.errorbar(Llist["rentails"],
            MassInf, MassErr, marker=marker, label=r"$E_T=\infty$")

    sigma = MassErr
    # sigma = None
    popt, pcov = curve_fit(Massfun, Llist["rentails"], MassInf.ravel(),
            bounds=boundsMass, method=method, sigma=sigma)
    xdata = scipy.linspace(xmin, xmax, 100)
    ax.plot(xdata, Massfun(xdata, *popt))
    msg = [
            r"$m_{{ph}} = {:.7f} \pm {:.7f}$".format(popt[0],np.sqrt(pcov[0,0])),
            r"$b = {:.7f} \pm {:.7f}$".format(popt[1],np.sqrt(pcov[1,1])),
            r"$c = {:.7f} \pm {:.7f}$".format(popt[2],np.sqrt(pcov[2,2]))
            ]
    for i, m in enumerate(msg):
        ax.text(0.8, 0.2-i*0.05, msg[i], horizontalalignment='center',
            verticalalignment='center', fontsize=13, transform=ax.transAxes)



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

# VACUUM ENERGY DENSITY
title = r"$g$={:.1f}, \quad $\alpha$={}, \quad f(x)={}".format(g, alpha,fvacStr)
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$E_0/L$")
ymargin = 10**(-5)
plt.xlim(xmin, xmax)
plt.ylim(ymin[1]-ymargin, ymax[1]+ymargin)
plt.legend(loc=2)
plt.savefig("extrLambdavsL_"+fname)

# MASS
title = r"$g$={:.1f}, \quad $\alpha$={}, \quad f(x)={}".format(g, alpha,fmassStr)
plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$L$")
plt.ylabel(r"$M$")
plt.xlim(xmin, xmax)
plt.ylim(ymin[-1]-ymargin, ymax[-1]+ymargin)
plt.legend(loc=1)
plt.savefig("extrMassvsL_"+fname)
