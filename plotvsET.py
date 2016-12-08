import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
from cycler import cycler
import database
from sys import exit

output = "png"
renlist = ("raw", "rentails", "renloc")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

klist = (1,-1)

neigs = 1


# Ratio between EL and ET
ratioELET = 1.5
# Ratio between ELp and ET
ratioELpET = 1.5
# Ratio between ELpp and ELp
ratioELppELp = 1.5

maxntails = 300


def fignum(k):
    if k==1:
        return 1
    elif k==-1:
        return 2


def plotvsET(ETlist):

    xlist = ETlist

    db = database.Database()

    spectrum = {k:{ren:[] for ren in renlist} for k in klist}
    mass = {ren:[] for ren in renlist}
    ntails = {k:[] for k in klist}

    for k in klist:

        for ET in ETlist:

            for ren in renlist:
                # print("ET=",ET)

                EL = ratioELET*ET
                ELp = ratioELpET*ET
                ELpp = ratioELppELp*ELp

                approxQuery = {"g":g, "L":L, "ET":ET}
                exactQuery = {"k": k, "ren":ren}
                boundQuery = {}

                if ren=="rentails":
                    exactQuery["maxntails"] = maxntails
                    exactQuery["tailsComputedAtET"] = ETmax
                    approxQuery["EL"] = EL
                    approxQuery["ELp"] = ELp
                    approxQuery["ELpp"] = ELpp


                try:
                    spectrum[k][ren].append(db.getObjList('spec', exactQuery, approxQuery,
                        boundQuery)[0])

                    if ren=="rentails":
                        ntails[k].append(db.getObjList('ntails', exactQuery, approxQuery,
                            boundQuery)[0])

                except IndexError:
                    print("Not found:", exactQuery, approxQuery)
                    exit(-1)


    # Convert to array
    for k in klist:
        for ren in renlist:
            spectrum[k][ren] = array(spectrum[k][ren])

    # Mass
    if -1 in klist and 1 in klist:
        for ren in renlist:
            mass[ren] = spectrum[-1][ren][:,0]-spectrum[1][ren][:,0]

    # SPECTRUM
    for k in klist:
        plt.figure(fignum(k))
        for ren in renlist:
            for i in range(neigs):
                data = spectrum[k][ren][:,i]
                plt.plot(xlist, data, label="ren="+ren)

    # MASS
    if -1 in klist and 1 in klist:
        plt.figure(3)
        for ren in renlist:
            data = mass[ren]
            plt.plot(xlist, data, label="ren="+ren)


    # Number of tails
    plt.figure(4)
    for k in klist:
        y = array(ntails[k])
        plt.plot(xlist,y,label="k="+str(k))


argv = sys.argv


if len(argv) < 5:
    print(argv[0], "<L> <g> <ETmin> <ETmax>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ETmin = float(argv[3])
ETmax = float(argv[4])

ETlist = scipy.linspace(ETmin, ETmax, (ETmax-ETmin)*2+1)
print("ETlist:", ETlist)


params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))


plotvsET(ETlist)



title = r"$g$={0:.1f}, $L$={1:.1f}, maxntails={2},"\
            "$E_L/E_T$={3:.1f}, $E_L'/E_T$={4:.1f}, $E_L''/E_T={5:.1f}$"\
            .format(g,L,maxntails,ratioELET,ratioELpET,ratioELppELp)
fname = "g={0:.1f}_L={1:.1f}_maxntails={2}.{3}_"\
            "ELET={3}_ELpET={4}_ELppELp={5}.{6}"\
            .format(g,L,maxntails,ratioELET,ratioELpET,ratioELppELp,output)
loc = "upper right"


# Even eigenvalues
if 1 in klist:
    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_{T}$")
    plt.ylabel(r"$E_i$ even")
    plt.legend(loc=loc)
    plt.savefig("evenvsET_"+fname)


# Odd eigenvalues
if -1 in klist:
    plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_{T}$")
    plt.ylabel(r"$E_i$ odd")
    plt.legend(loc=loc)
    plt.savefig("oddvsET_"+fname)


# Mass
if 1 in klist and -1 in klist:
    plt.figure(3, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_{T}$")
    plt.ylabel(r"$m_{\rm ph}$")
    plt.legend(loc=loc)
    plt.savefig("massvsET_"+fname)


# Plot of the number of tails
title = r"$g$={0:.1f}, $L$={1:.1f}, maxntails={2}".format(g,L,maxntails)
fname = "g={0:.1f}_L={1:.1f}_maxntails={2}.{3}".format(g,L,maxntails,output)

plt.figure(4, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_{T}$")
plt.ylabel(r"numtails")
plt.legend(loc=loc)
plt.savefig("tailsvsET_"+fname)
