import sys
import matplotlib.pyplot as plt
import scipy
import math
import json
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
import database
import phi1234

output = "svg"

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def main(argv):
    args = " <g> <Emax> <fname>"
    if len(argv) < 4:
        print(argv[0], args)
        return -1

    g = float(argv[1])
    Emax = float(argv[2])
    fname = argv[3]

    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)


    a = phi1234.Phi1234()
    a.loadMatrix(fname)
    L = a.L
    a.buildBasis(k=1, Emax=Emax)

    db = database.Database()
    exactQuery = {"ren":"raw", "k":1}
    approxQuery = {"g":g, "Emax":Emax, "L":L}
    eigv = db.getObjList("eigv", exactQuery=exactQuery, approxQuery=approxQuery)[0]

    # Select only basis coefficients with 2 particles
    vacuumEigv = eigv[0]
    wf = [vacuumEigv[i]  for i in range(len(vacuumEigv)) if a.basis[1][i].occN()==2]
    print(wf)


    plt.figure(1)
    plt.plot(range(len(wf)), wf, "o")
            # linewidth=0.7,
            # dashes = [3,2], marker='.', markersize=3, color=evenColor, label=label("raw"))


    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    # #plt.xlim(min(xList)-0.01, max(xList)+0.01)
    plt.title("L={:.1f}, Emax={:.1f}, g={:.1f}".format(L,Emax,g))
    # plt.xlabel(r"$L$")
    # plt.ylabel(r"$\Lambda$")
    # plt.legend(loc='lower right', prop={'size':10})
    plt.savefig("wf_L={:.0f}_Emax={:.0f}_g={:.0f}.{:s}".format(L,Emax,g,output))


if __name__ == "__main__":
    main(sys.argv)
