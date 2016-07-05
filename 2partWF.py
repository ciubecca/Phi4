import sys
import matplotlib.pyplot as plt
import scipy
import math
import json
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
import database
import statefuncs

output = "pdf"

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def main(argv):
    args = " <g1,g2,...,gN> <Emax> <L>"
    if len(argv) < 4:
        print(argv[0], args)
        return -1

    gList = [float(x) for x in argv[1].split(",")]
    Emax = float(argv[2])
    L = float(argv[3])

    plt.figure(1)
    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)

    basis = statefuncs.Basis(m=1, L=L, Emax=Emax, k=1)

    for g in gList:
        db = database.Database()
        exactQuery = {"ren":"raw", "k":1}
        approxQuery = {"g":g, "Emax":Emax, "L":L}
        eigv = db.getObjList("eigv", exactQuery=exactQuery, approxQuery=approxQuery)[0]

        # Select only basis coefficients with 2 particles
        vacuumEigv = eigv[0]
        wf = array([vacuumEigv[i]  for i in range(len(vacuumEigv)) if basis[i].occN()==2])

        # Normalization: multiply 2 particles at rest by sqrt(2)
        # Takes into account normalization of states and Bose symmetry
        # XXX check the normalization
        wf[0] = wf[0]*sqrt(2)

        # Since the wave function will be O(g^2) in the perturbative limit, rescale by this parameter
        wf = wf/g**2

        klist = array(range(len(wf)))*(2*pi)/L
        plt.plot(klist,wf, label="g={:.1f}".format(g))
            # linewidth=0.7,
            # dashes = [3,2], marker='.', markersize=3, color=evenColor, label=label("raw"))


    # plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    # #plt.xlim(min(xList)-0.01, max(xList)+0.01)
    plt.title("L={:.1f}, Emax={:.1f}".format(L,Emax))
    plt.xlabel(r"$k$")
    # plt.ylim(min(wf)-0.01,max(wf)+0.01)
    # plt.ylabel(r"$\Lambda$")
    plt.legend(loc='lower right', prop={'size':10})
    plt.savefig("wf_L={:.0f}_Emax={:.0f}.{:s}".format(L,Emax,output))


if __name__ == "__main__":
    main(sys.argv)
