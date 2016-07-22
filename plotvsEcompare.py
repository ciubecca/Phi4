import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
import database
import finiteVolH

output = "pdf"
renlist = ("raw", "renlocal")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def main(argv):
    args = "<L> <g>"
    if len(argv) < 3:
        print("{0} {1}".format(argv[0],args))
        return -1

    L = float(argv[1])
    g = float(argv[2])
    occmax = 4
    Emaxbar = 30

    plotWithTails(L,g,occmax,Emaxbar)
    plotWithoutTails(L,g,occmax)

    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    #plt.xlim(min(xList)-0.01, max(xList)+0.01)
    plt.title(r"$g$={0:.2f}, $L$={1:.2f}, $\bar{{E}}_{{\rm max}}$={2:.2f}, $n_{{\rm max}}$={3:d}".format(g,L,Emaxbar,occmax))
    plt.xlabel(r"$E_{{\rm max}}$")
    plt.ylabel(r"$E_0$")
    plt.legend()
    plt.savefig("fig_E0vsEcompare_g={0:.2f}_L={1:.2f}_nmax={2:d}.{3}".format(g,L,occmax,output))



def plotWithTails(L, g, occmax, Emaxbar):
    # Hardcoded parameters
    neigs = 1
    Elist = scipy.linspace(5, 28, 24)
    print(Elist)

    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)

    db = database.Database()

    exactQuery = {"k":1, "occmax":occmax}
    approxQuery = {"g":g, "L":L, "Emaxbar":Emaxbar}

    E0 = {}
    for ren in renlist:
        E0[ren] = []

        for Emax in Elist:
            exactQuery["ren"] = ren
            approxQuery["Emax"] = Emax
            E0[ren].append(db.getObjList('spec', exactQuery, approxQuery)[0][0])


    # VACUUM ENERGY
    plt.figure(1)

    data = E0["raw"]
    plt.plot(Elist, data, linewidth=1., color="b",
            dashes = [3,2], label="raw")

    data = E0["renlocal"]
    plt.plot(Elist, data, linewidth=1., color="r",
            dashes = [3,2], label="renlocal")



def plotWithoutTails(L, g, occmax):
    # Hardcoded parameters
    neigs = 1
    Elist = scipy.linspace(5, 28, 24)
    print(Elist)

    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)

    db = database.Database()

    exactQuery = {"k":1, "occmax":occmax}
    approxQuery = {"g":g, "L":L}

    E0 = {}
    for ren in renlist:
        E0[ren] = []

        for Emax in Elist:
            exactQuery["ren"] = ren
            approxQuery["Emax"] = Emax
            approxQuery["Emaxbar"] = Emax
            E0[ren].append(db.getObjList('spec', exactQuery, approxQuery)[0][0])


    # VACUUM ENERGY
    plt.figure(1)

    data = E0["raw"]
    plt.plot(Elist, data, linewidth=1., color="b",
            dashes = [4,1], marker='+', markersize=3, label="raw no tails")

    data = E0["renlocal"]
    plt.plot(Elist, data, linewidth=1., color="r",
            dashes = [4,1], marker='+', markersize=3, label="renlocal no tails")


if __name__ == "__main__":
    main(sys.argv)
