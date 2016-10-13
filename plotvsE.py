# This files generates plots of Vacuum and mass eigenvalues from the database
# It should be called as:
# plotvsE L g ETmin ETmax
# For instance:
# plotvsE 10 1 10 20
# Plots all the points for L=10, g=1, and ET = [10, 10.5, 11, 11.5, ..., 20]


import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
import database

output = "pdf"
renlist = ("raw", "renloc", "rentails")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

neigs = 10

def EL(ET):
    return ET*2.5

def plotvsE(Elist):

    db = database.Database()

    exactQuery = {}
    approxQuery = {"g":g, "L":L}

    E0 = {}
    E1 = {}
    for ren in renlist:
        E0[ren] = []
        E1[ren] = []

        for ET in ETlist:
            exactQuery["ren"] = ren
            approxQuery["ET"] = ET
            exactQuery["k"] = 1
            if ren=="renloc":
                approxQuery["EL"] = ET
            elif ren=="rentails":
                approxQuery["EL"] = EL(ET)
                # approxQuery["EL"] = max(ETlist)+ELETdiff
            E0[ren].append(db.getObjList('spec', exactQuery, approxQuery)[0][0])

            exactQuery["k"] = -1
            E1[ren].append(db.getObjList('spec', exactQuery, approxQuery)[0][0])


    # VACUUM ENERGY
    plt.figure(1)

    linewidth = 1
    dashes = [4,4]
    marker = "."
    markersize = 6

    data = E0["raw"]
    plt.plot(Elist, data, linewidth=linewidth, color="b", marker=marker,
            markersize=markersize, dashes = dashes, label="raw")

    data = E0["renloc"]
    plt.plot(Elist, data, linewidth=linewidth, color="r", marker=marker,
            markersize=markersize, dashes = dashes, label="renloc")

    data = E0["rentails"]
    plt.plot(Elist, data, linewidth=linewidth, color="k", marker=marker,
            markersize=markersize, dashes = dashes, label="rentails")

    # MASS
    plt.figure(2)

    data = array(E1["raw"])-array(E0["raw"])
    plt.plot(Elist, data, linewidth=linewidth, color="b", marker=marker,
            markersize=markersize, dashes = dashes, label="raw")

    data = array(E1["renloc"])-array(E0["renloc"])
    plt.plot(Elist, data, linewidth=linewidth, color="r", marker=marker,
            markersize=markersize, dashes = dashes, label="renloc")

    data = array(E1["rentails"])-array(E0["rentails"])
    plt.plot(Elist, data, linewidth=linewidth, color="k", marker=marker,
            markersize=markersize, dashes = dashes, label="rentails")




argv = sys.argv
if len(argv) < 5:
    print(argv[0], "<L> <g> <ETmin> <ETmax>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ETmin = float(argv[3])
ETmax = float(argv[4])
# ELETdiff = float(argv[5])

ETlist = scipy.linspace(ETmin, ETmax, 2*(ETmax-ETmin)+1)
print("ETlist:", ETlist)


params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plotvsE(ETlist)

plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
#plt.xlim(min(xList)-0.01, max(xList)+0.01)
plt.title(r"$g$={0:.1f}, $L$={1:.1f}, $E_L=2.5 E_T$".format(g,L))
# plt.title(r"$g$={0:.1f}, $L$={1:.1f}, $E_L$={2:.1f}".format(g,L,EL))
plt.xlabel(r"$E_T$")
plt.ylabel(r"$E_0$")
plt.legend(loc="upper right")


plt.savefig("figs/fig_E0vsET_g={0:.1f}_L={1:.1f}_EL=2.5*ET.{2}"
        .format(g,L,output))
# plt.savefig("figs/fig_E0vsET_g={0:.1f}_L={1:.1f}_EL={2:.1f}.{3}"
        # .format(g,L,EL,output))

plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
#plt.xlim(min(xList)-0.01, max(xList)+0.01)
plt.title(r"$g$={0:.1f}, $L$={1:.1f}, $E_L=2.5 E_T$".format(g,L))
# plt.title(r"$g$={0:.1f}, $L$={1:.1f}, $E_L$={2:.1f}".format(g,L,EL))
plt.xlabel(r"$E_T$")
plt.ylabel(r"$E_1-E_0$")
plt.legend(loc="upper right")

plt.savefig("figs/fig_MvsET_g={0:.1f}_L={1:.1f}_EL=2.5*ET.{2}"
        .format(g,L,output))
# plt.savefig("figs/fig_MvsET_g={0:.1f}_L={1:.1f}_EL={2:.1f}.{3}"
        # .format(g,L,EL,output))


