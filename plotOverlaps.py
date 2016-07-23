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
    args = "<Emax> <g>"
    if len(argv) < 3:
        print("{0} {1}".format(argv[0],args))
        return -1

    Emax = float(argv[1])
    g = float(argv[2])

    # Hardcoded parameters
    neigs = 1
    Llist = scipy.linspace(5, 15, 21)
    print(Llist)

    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)

    db = database.Database()

    exactQuery = {"k":1}
    approxQuery = {"g":g, "Emax":Emax}

    E0 = {}
    v0 = {}
    for ren in renlist:
        E0[ren] = []
        v0[ren] = []

        for L in Llist:
            exactQuery["ren"] = ren
            approxQuery["L"] = L
            E0[ren].append(db.getObjList('spec', exactQuery, approxQuery)[0][0])
            v0[ren].append(db.getObjList('eigv', exactQuery, approxQuery)[0][0])


    # VACUUM ENERGY
    plt.figure(1)

    data = E0["renlocal"]/Llist
    plt.plot(Llist, data, linewidth=1.,
            dashes = [4,1], marker='+', markersize=3)


    # OVERLAP WITH PERTURBATIVE VACUUM
    plt.figure(2)
    data = [x[0]**2. for x in v0["renlocal"]]
    plt.plot(Llist, data, linewidth=0.7,
            dashes = [3,2], marker='.', markersize=3)

#############################################################################

    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    #plt.xlim(min(xList)-0.01, max(xList)+0.01)
    plt.title(r"g={0:.2f}, $E_{{\rm max}}$={1:.2f}".format(g,Emax))
    plt.xlabel(r"$L$")
    plt.ylabel(r"$E_0/L$")
    plt.savefig("fig_E0vsL_g={0:.2f}_Emax={1:.2f}.{2}".format(g,Emax,output))

    plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    #plt.xlim(min(xList)-0.01, max(xList)+0.01)
    plt.title(r"g={0:.2f}, $E_{{\rm max}}$={1:.2f}".format(g,Emax))
    plt.xlabel(r"$L$")
    plt.ylabel(r"$\langle 0 \mid \psi_0 \rangle$")
    plt.savefig("fig_Psi0vsL_g={0:.2f}_Emax={1:.2f}.{2}".format(g,Emax,output))


if __name__ == "__main__":
    main(sys.argv)
