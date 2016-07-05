import sys
import matplotlib.pyplot as plt
import scipy
import math
import json
from scipy import pi, log, array, sqrt
from math import factorial
from matplotlib import rc
import database
import statefuncs

output = "pdf"

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


""" Takes the a coefficient of the discretized wafe functions and renormalizes it
    according to the basis element.
    Works for general particle number """
# XXX check this
def normalizeWF(c, v):
    # Bose symmetry
    c *= scipy.prod([factorial(n) for n in v])/factorial(v.occN())

    # Normalization of Fock space states
    c *= 1/sqrt(scipy.prod([factorial(n) for n in v]))

    # Spatial parity symmetry
    if v.isParityEigenstate() == False:
        c *= 1/sqrt(2)

    return c

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

    basis = statefuncs.Basis(m=1, L=L, Emax=Emax, k=-1)
    # Select only 3 particles basis states
    indexList = [i for i in range(len(basis)) if basis[i].occN()==3]
    basis3p = [basis[i] for i in indexList]

    ydata = []

    for g in gList:
        db = database.Database()
        exactQuery = {"ren":"raw", "k":-1}
        approxQuery = {"g":g, "Emax":Emax, "L":L}
        eigv = db.getObjList("eigv", exactQuery=exactQuery, approxQuery=approxQuery)[0]

        # Select only coefficients of 3 particles basis states
        wf = [eigv[0][i] for i in indexList]

        # Renormalize wave function
        wf = array([normalizeWF(c,v) for c,v in zip(wf, basis3p)])

        print(wf)

        # Since the wave function will be O(g^2) in the perturbative limit, rescale by this parameter
        # wf = wf/g**2

        # Change sign to compare wave functions
        # if wf[0]<0:
            # wf = -wf

        # klist = array(range(len(wf)))*(2*pi)/L
        # plt.plot(klist,wf, label="g={:.3f}".format(g))
            # # linewidth=0.7,
            # # dashes = [3,2], marker='.', markersize=3, color=evenColor, label=label("raw"))
        # ydata.extend(wf)

    # plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    # #plt.xlim(min(xList)-0.01, max(xList)+0.01)
    # plt.title("L={:.1f}, Emax={:.1f}".format(L,Emax))
    # plt.xlabel(r"$k$")
    # plt.ylabel(r"$\psi(k)/g^2$")
    # plt.ylim(min(ydata)-0.01,max(ydata)+0.01)
    # plt.xlim(min(klist),max(klist))
    # # plt.ylabel(r"$\Lambda$")
    # plt.legend(loc='upper right', prop={'size':10})
    # plt.savefig("wf_L={:.0f}_Emax={:.0f}.{:s}".format(L,Emax,output))


if __name__ == "__main__":
    main(sys.argv)
