# This files generates plots of Vacuum and mass eigenvalues from the database
# It should be called as:
# plotvsE L g ETmin ETmax
# For instance:
# plotvsE 10 1 10 20
# Plots all the points for L=10, g=1, and ET = [10, 10.5, 11, 11.5, ..., 20]

minoverlaplist = [10**(-2)]

import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
from cycler import cycler
import database

output = "png"
# renlist = ("raw", "renloc", "rentails")
renlist = ("rentails",)

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

neigs = 10


def plotvsEL2(minoverlap, EL2list):

    xlist = EL2list

    db = database.Database()

    exactQuery = {}
    approxQuery = {"g":g, "L":L}

    E0 = {}
    E1 = {}
    M = {}
    for ren in renlist:
        E0[ren] = []
        E1[ren] = []

        for EL2 in EL2list:
            exactQuery["ren"] = ren
            approxQuery["ET"] = ET
            approxQuery["minoverlap"] = minoverlap
            approxQuery["EL"] = EL2
            approxQuery["EL3"] = min(EL2list)

            exactQuery["k"] = 1
            E0[ren].append(db.getObjList('spec', exactQuery, approxQuery)[0][0])

            exactQuery["k"] = -1
            E1[ren].append(db.getObjList('spec', exactQuery, approxQuery)[0][0])


    label = "minoverlap={}".format(minoverlap)

    # VACUUM ENERGY
    plt.figure(1)

    # data = E0["raw"]
    # plt.plot(Elist, data, linewidth=linewidth, color="b", marker=marker,
            # markersize=markersize, dashes = dashes, label="raw")

    # data = E0["renloc"]
    # plt.plot(Elist, data, linewidth=linewidth, color="r", marker=marker,
            # markersize=markersize, dashes = dashes, label="renloc")


    data = E0["rentails"]
    plt.plot(xlist, data, label=label)



    # MASS
    plt.figure(2)

    # data = array(E1["raw"])-array(E0["raw"])
    # plt.plot(Elist, data, linewidth=linewidth, color="b", marker=marker,
            # markersize=markersize, dashes = dashes, label="raw")

    # data = array(E1["renloc"])-array(E0["renloc"])
    # plt.plot(Elist, data, linewidth=linewidth, color="r", marker=marker,
            # markersize=markersize, dashes = dashes, label="renloc")

    data = array(E1["rentails"])-array(E0["rentails"])
    plt.plot(xlist, data, label=label)



    plt.figure(3)
    data = array(E1["rentails"])
    plt.plot(xlist, data, label=label)



argv = sys.argv


if len(argv) < 6:
    print(argv[0], "<L> <g> <ET> <EL2min> <EL2max>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
EL2min = float(argv[4])
EL2max = float(argv[5])



EL2list = scipy.linspace(EL2min, EL2max, (EL2max-EL2min)+1)
print("EL2list:", EL2list)


params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))


for minoverlap in minoverlaplist:
    plotvsEL2(minoverlap, EL2list)


title = r"$g$={0:.1f}, $L$={1:.1f}, $E_T$={2:.1f}".format(g,L,ET)
fname = "noLoc_g={0:.1f}_L={1:.1f}_ET={2:.1f}.{3}".format(g,L,ET,output)
loc = "lower right"

plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
#plt.xlim(min(xList)-0.01, max(xList)+0.01)
plt.title(title)
# plt.title(r"$g$={0:.1f}, $L$={1:.1f}, $E_L$={2:.1f}".format(g,L,EL))
plt.xlabel(r"$E_{L 2}$")
plt.ylabel(r"$E_0$")
plt.legend(loc=loc)


plt.savefig("figs/fig_E0vsEL2_"+fname)
# plt.savefig("figs/fig_E0vsET_g={0:.1f}_L={1:.1f}_EL={2:.1f}.{3}"
        # .format(g,L,EL,output))

plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
#plt.xlim(min(xList)-0.01, max(xList)+0.01)
plt.title(title)
# plt.title(r"$g$={0:.1f}, $L$={1:.1f}, $E_L$={2:.1f}".format(g,L,EL))
plt.xlabel(r"$E_{L 2}$")
plt.ylabel(r"$E_1-E_0$")
plt.legend(loc=loc)

plt.savefig("figs/fig_MvsEL2_"+fname)
# plt.savefig("figs/fig_MvsET_g={0:.1f}_L={1:.1f}_EL={2:.1f}.{3}"
        # .format(g,L,EL,output))



plt.figure(3, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
#plt.xlim(min(xList)-0.01, max(xList)+0.01)
plt.title(title)
# plt.title(r"$g$={0:.1f}, $L$={1:.1f}, $E_L$={2:.1f}".format(g,L,EL))
plt.xlabel(r"$E_{L 2}$")
plt.ylabel(r"$E_1$")
plt.legend(loc=loc)

<<<<<<< HEAD
=======

>>>>>>> numRenNoloc3
plt.savefig("figs/fig_E1vsEL2_"+fname)
# plt.savefig("figs/fig_MvsET_g={0:.1f}_L={1:.1f}_EL={2:.1f}.{3}"
        # .format(g,L,EL,output))