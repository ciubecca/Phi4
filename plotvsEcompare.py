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

k = 1
neigs = 1
Emin = 5
Emaxbar = None
occmax = None
L = None
g = None

def main(argv):
    args = "<L> <Emaxbar> <g> <occmax>"
    if len(argv) < 4:
        print("{0} {1}".format(argv[0],args))
        return -1

    global L
    global Emaxbar
    global g
    global occmax

    L = float(argv[1])
    Emaxbar = float(argv[2])
    g = float(argv[3])
    try:
        occmax = int(argv[4])
    except IndexError:
        occmax = None

    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)

    Elist = scipy.linspace(Emin, Emaxbar-1, Emaxbar-Emin)
    print(Elist)

    plotvsE(Elist, tails=True)
    plotvsE(Elist, tails=False)

    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    #plt.xlim(min(xList)-0.01, max(xList)+0.01)
    plt.title(
      r"$g$={0:.2f}, $L$={1:.2f}, $\bar{{E}}_{{\rm max}}$={2:.2f},$n_{{\rm max}}$={3:d}"
        .format(g,L,Emaxbar,occmax))
    plt.xlabel(r"$E_{{\rm max}}$")
    plt.ylabel(r"$E_0$")
    plt.legend(loc="lower right")
    plt.savefig("figs/fig_E0vsEcompare_g={0:.2f}_L={1:.2f}_nmax={2:d}.{3}"
            .format(g,L,occmax,output))



def plotvsE(Elist, tails):

    db = database.Database()

    exactQuery = {"k":k, "occmax":occmax}
    approxQuery = {"g":g, "L":L}

    E0 = {}
    for ren in renlist:
        E0[ren] = []

        for Emax in Elist:
            if tails==True:
                cutoff = Emaxbar
            else:
                cutoff = Emax

            exactQuery["ren"] = ren
            approxQuery["Emax"] = Emax
            approxQuery["Emaxbar"] = cutoff
            E0[ren].append(db.getObjList('spec', exactQuery, approxQuery)[0][0])


    # VACUUM ENERGY
    plt.figure(1)

    if tails==True:
        s = ""
        linewidth = .8
        dashes = [4,1]
        marker = '+'
        markersize = 4
    else:
        s = "o"
        linewidth = 1
        dashes = [4,4]
        marker = None
        markersize = None

    data = E0["raw"]
    plt.plot(Elist, data, linewidth=linewidth, color="b", marker=marker,
            markersize=markersize, dashes = dashes, label="raw w/"+s+" tails")

    if tails:
        plt.axhline(y=data[-1], color='k')
        # plt.axhline(y=data[-1], xmin=min(Elist), xmax=max(Elist), linewidth=2, color = 'k')

    data = E0["renlocal"]
    plt.plot(Elist, data, linewidth=linewidth, color="r", marker=marker,
            markersize=markersize, dashes = dashes, label="ren w/"+s+" tails")




if __name__ == "__main__":
    main(sys.argv)
